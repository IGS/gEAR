#!/opt/bin/python3

"""
stats.cgi — Google Analytics 4 statistics endpoint

Query Parameters
----------------
days : int, optional
    Number of past days to include in the report window.
    Default: 30. Minimum: 1. Maximum: 365.
    Example: ?days=7

start_date : str, optional
    Absolute start date in YYYY-MM-DD format.
    If provided, it overrides the start date implied by `days`.
    Example: ?start_date=2026-01-01

end_date : str, optional
    Absolute end date in YYYY-MM-DD format.
    If provided, it overrides the default end date of "today".
    Example: ?end_date=2026-01-31

top_n : int, optional
    Maximum number of rows to return for ranked reports
    (top pages, traffic sources, countries, top events).
    Default: 10. Minimum: 1. Maximum: 50.
    Example: ?top_n=5

realtime : str, optional
    Whether to include real-time active-user data.
    Truthy values: any value other than "0", "false", or "no".
    Default: 1 (enabled).
    Example: ?realtime=0

Response Shape (JSON)
---------------------
{
  "ok": true,
  "range": {
    "days": <int>,
        "startDate": "<effective start date>",
        "endDate": "<effective end date>"
  },
  "summary": {
    "activeUsers": <int>,
    "sessions": <int>,
    "screenPageViews": <int>,
    "engagedSessions": <int>,
    "eventCount": <int>
  },
  "timeseries": [
    {"date": "YYYYMMDD", "activeUsers": <int>, "sessions": <int>, "screenPageViews": <int>},
    ...
  ],
  "topPages": [
    {"pagePath": <str>, "screenPageViews": <int>, "activeUsers": <int>, "sessions": <int>},
    ...
  ],
  "trafficSources": [
    {"sessionSourceMedium": <str>, "sessions": <int>, "activeUsers": <int>},
    ...
  ],
  "countries": [
    {"country": <str>, "activeUsers": <int>, "sessions": <int>, "screenPageViews": <int>},
    ...
  ],
  "topEvents": [
    {"eventName": <str>, "eventCount": <int>, "activeUsers": <int>},
    ...
  ],
  "realtime": {                  // only present when realtime=1 (default)
    "activeUsersNow": <int>,
    "byCountry": [
      {"country": <str>, "activeUsers": <int>},
      ...
    ]
  }
}
"""

import json
import os
import re
import sys
import traceback
from datetime import datetime
from urllib.parse import parse_qs

# from gear.serverconfig import ServerConfig
from google.analytics.data_v1beta import BetaAnalyticsDataClient
from google.analytics.data_v1beta.types import (
    DateRange,
    Dimension,
    Metric,
    OrderBy,
    RunRealtimeReportRequest,
    RunReportRequest,
)

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
# server_config = ServerConfig().parse()
# GA4_PROPERTY_ID = os.environ.get("GA4_PROPERTY_ID", "").strip()
GA4_PROPERTY_ID = "312199300"

# Optional: restrict CORS if this CGI will be called from JS on another origin.
# Example: "https://umgear.org"
ALLOWED_ORIGIN = os.environ.get("GEAR_STATS_ALLOWED_ORIGIN", "").strip()

DEFAULT_DAYS = 30
MAX_DAYS = 365

# -----------------------------------------------------------------------------
# Utility helpers
# -----------------------------------------------------------------------------


def print_json_response(payload, status="200 OK"):
    print(f"Status: {status}")
    print("Content-Type: application/json")
    print("Cache-Control: no-store, no-cache, must-revalidate, max-age=0")
    print("Pragma: no-cache")
    if ALLOWED_ORIGIN:
        print(f"Access-Control-Allow-Origin: {ALLOWED_ORIGIN}")
    print()
    print(json.dumps(payload, indent=2))


def error_response(message, status="500 Internal Server Error", details=None):
    payload = {
        "ok": False,
        "error": message,
    }
    if details is not None:
        payload["details"] = details
    print_json_response(payload, status=status)
    sys.exit(0)


def parse_int(value, default, minimum=None, maximum=None):
    try:
        ivalue = int(value)
    except (TypeError, ValueError):
        return default

    if minimum is not None and ivalue < minimum:
        ivalue = minimum
    if maximum is not None and ivalue > maximum:
        ivalue = maximum
    return ivalue


def safe_int(value):
    try:
        return int(value)
    except (TypeError, ValueError):
        try:
            return int(float(value))
        except (TypeError, ValueError):
            return 0


def parse_iso_date_param(value, name):
    if value is None:
        return None

    value = value.strip()
    if not value:
        return None

    if not re.fullmatch(r"\d{4}-\d{2}-\d{2}", value):
        error_response(
            f"Invalid {name} format. Expected YYYY-MM-DD.",
            status="400 Bad Request",
        )

    try:
        datetime.strptime(value, "%Y-%m-%d")
    except ValueError:
        error_response(
            f"Invalid {name} value. Expected a real calendar date in YYYY-MM-DD format.",
            status="400 Bad Request",
        )

    return value


def get_query_params():
    raw = os.environ.get("QUERY_STRING", "")
    qs = parse_qs(raw, keep_blank_values=False)

    days = parse_int(
        qs.get("days", [str(DEFAULT_DAYS)])[0],
        default=DEFAULT_DAYS,
        minimum=1,
        maximum=MAX_DAYS,
    )

    top_n = parse_int(
        qs.get("top_n", ["10"])[0],
        default=10,
        minimum=1,
        maximum=50,
    )

    include_realtime = qs.get("realtime", ["1"])[0].lower() not in ("0", "false", "no")
    start_date = parse_iso_date_param(qs.get("start_date", [None])[0], "start_date")
    end_date = parse_iso_date_param(qs.get("end_date", [None])[0], "end_date")

    if start_date and end_date:
        start_dt = datetime.strptime(start_date, "%Y-%m-%d")
        end_dt = datetime.strptime(end_date, "%Y-%m-%d")
        if start_dt > end_dt:
            error_response(
                "Invalid date range. start_date must be less than or equal to end_date.",
                status="400 Bad Request",
            )

    return {
        "days": days,
        "top_n": top_n,
        "include_realtime": include_realtime,
        "start_date": start_date,
        "end_date": end_date,
    }


def build_client():
    return BetaAnalyticsDataClient()


def run_report(
    client,
    property_id,
    dimensions,
    metrics,
    start_date,
    end_date,
    limit=10,
    order_by_metric=None,
):
    request = RunReportRequest(
        property=f"properties/{property_id}",
        dimensions=[Dimension(name=d) for d in dimensions],
        metrics=[Metric(name=m) for m in metrics],
        date_ranges=[DateRange(start_date=start_date, end_date=end_date)],
        limit=limit,
        order_bys=(
            [
                OrderBy(
                    metric=OrderBy.MetricOrderBy(metric_name=order_by_metric), desc=True
                )
            ]
            if order_by_metric
            else []
        ),
    )
    return client.run_report(request=request)


def run_realtime_report(client, property_id, dimensions, metrics, limit=10):
    request = RunRealtimeReportRequest(
        property=f"properties/{property_id}",
        dimensions=[Dimension(name=d) for d in dimensions],
        metrics=[Metric(name=m) for m in metrics],
        limit=limit,
    )
    return client.run_realtime_report(request=request)


# -----------------------------------------------------------------------------
# Report builders
# -----------------------------------------------------------------------------


def get_summary(client, property_id, start_date, end_date):
    response = run_report(
        client=client,
        property_id=property_id,
        dimensions=[],
        metrics=[
            "activeUsers",
            "sessions",
            "screenPageViews",
            "engagedSessions",
            "eventCount",
        ],
        start_date=start_date,
        end_date=end_date,
        limit=1,
    )

    if not response.rows:
        return {
            "activeUsers": 0,
            "sessions": 0,
            "screenPageViews": 0,
            "engagedSessions": 0,
            "eventCount": 0,
        }

    row = response.rows[0]
    return {
        "activeUsers": safe_int(row.metric_values[0].value),
        "sessions": safe_int(row.metric_values[1].value),
        "screenPageViews": safe_int(row.metric_values[2].value),
        "engagedSessions": safe_int(row.metric_values[3].value),
        "eventCount": safe_int(row.metric_values[4].value),
    }


def get_daily_timeseries(client, property_id, start_date, end_date, limit):
    response = run_report(
        client=client,
        property_id=property_id,
        dimensions=["date"],
        metrics=["activeUsers", "sessions", "screenPageViews"],
        start_date=start_date,
        end_date=end_date,
        limit=limit,
        order_by_metric=None,
    )

    rows = []
    for row in response.rows:
        rows.append(
            {
                "date": row.dimension_values[0].value,  # YYYYMMDD
                "activeUsers": safe_int(row.metric_values[0].value),
                "sessions": safe_int(row.metric_values[1].value),
                "screenPageViews": safe_int(row.metric_values[2].value),
            }
        )
    return rows


def get_top_pages(client, property_id, start_date, end_date, top_n):
    response = run_report(
        client=client,
        property_id=property_id,
        dimensions=["pagePathAndQueryString"],
        metrics=["screenPageViews", "activeUsers", "sessions"],
        start_date=start_date,
        end_date=end_date,
        limit=top_n,
        order_by_metric="screenPageViews",
    )

    rows = []
    for row in response.rows:
        rows.append(
            {
                "pagePath": row.dimension_values[0].value,
                "screenPageViews": safe_int(row.metric_values[0].value),
                "activeUsers": safe_int(row.metric_values[1].value),
                "sessions": safe_int(row.metric_values[2].value),
            }
        )
    return rows


def get_traffic_sources(client, property_id, start_date, end_date, top_n):
    response = run_report(
        client=client,
        property_id=property_id,
        dimensions=["sessionSourceMedium"],
        metrics=["sessions", "activeUsers"],
        start_date=start_date,
        end_date=end_date,
        limit=top_n,
        order_by_metric="sessions",
    )

    rows = []
    for row in response.rows:
        rows.append(
            {
                "sessionSourceMedium": row.dimension_values[0].value,
                "sessions": safe_int(row.metric_values[0].value),
                "activeUsers": safe_int(row.metric_values[1].value),
            }
        )
    return rows


def get_countries(client, property_id, start_date, end_date, top_n):
    response = run_report(
        client=client,
        property_id=property_id,
        dimensions=["country"],
        metrics=["activeUsers", "sessions", "screenPageViews"],
        start_date=start_date,
        end_date=end_date,
        limit=top_n,
        order_by_metric="activeUsers",
    )

    rows = []
    for row in response.rows:
        rows.append(
            {
                "country": row.dimension_values[0].value,
                "activeUsers": safe_int(row.metric_values[0].value),
                "sessions": safe_int(row.metric_values[1].value),
                "screenPageViews": safe_int(row.metric_values[2].value),
            }
        )
    return rows


def get_top_events(client, property_id, start_date, end_date, top_n):
    response = run_report(
        client=client,
        property_id=property_id,
        dimensions=["eventName"],
        metrics=["eventCount", "activeUsers"],
        start_date=start_date,
        end_date=end_date,
        limit=top_n,
        order_by_metric="eventCount",
    )

    rows = []
    for row in response.rows:
        rows.append(
            {
                "eventName": row.dimension_values[0].value,
                "eventCount": safe_int(row.metric_values[0].value),
                "activeUsers": safe_int(row.metric_values[1].value),
            }
        )
    return rows


def get_realtime(client, property_id, top_n):
    response = run_realtime_report(
        client=client,
        property_id=property_id,
        dimensions=["country"],
        metrics=["activeUsers"],
        limit=top_n,
    )

    by_country = []
    total_active_users = 0

    for row in response.rows:
        users = safe_int(row.metric_values[0].value)
        total_active_users += users
        by_country.append(
            {
                "country": row.dimension_values[0].value,
                "activeUsers": users,
            }
        )

    return {
        "activeUsersNow": total_active_users,
        "byCountry": by_country,
    }


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------


def main():
    if not GA4_PROPERTY_ID:
        error_response(
            "Missing GA4_PROPERTY_ID environment variable",
            status="500 Internal Server Error",
        )

    params = get_query_params()
    start_date = params["start_date"] or f"{params['days']}daysAgo"
    end_date = params["end_date"] or "today"
    timeseries_limit = params["days"] + 5

    if params["start_date"] or params["end_date"]:
        # Use a large cap when an explicit date parameter is used so daily rows
        # are not truncated by the default days-based limit.
        timeseries_limit = 10000

    if params["start_date"] and params["end_date"]:
        start_dt = datetime.strptime(params["start_date"], "%Y-%m-%d")
        end_dt = datetime.strptime(params["end_date"], "%Y-%m-%d")
        timeseries_limit = (end_dt - start_dt).days + 6

    try:
        client = build_client()

        payload = {
            "ok": True,
            "range": {
                "days": params["days"],
                "startDate": start_date,
                "endDate": end_date,
            },
            "summary": get_summary(client, GA4_PROPERTY_ID, start_date, end_date),
            "timeseries": get_daily_timeseries(
                client,
                GA4_PROPERTY_ID,
                start_date,
                end_date,
                timeseries_limit,
            ),
            "topPages": get_top_pages(
                client, GA4_PROPERTY_ID, start_date, end_date, params["top_n"]
            ),
            "trafficSources": get_traffic_sources(
                client, GA4_PROPERTY_ID, start_date, end_date, params["top_n"]
            ),
            "countries": get_countries(
                client, GA4_PROPERTY_ID, start_date, end_date, params["top_n"]
            ),
            "topEvents": get_top_events(
                client, GA4_PROPERTY_ID, start_date, end_date, params["top_n"]
            ),
        }

        if params["include_realtime"]:
            payload["realtime"] = get_realtime(client, GA4_PROPERTY_ID, params["top_n"])

        print_json_response(payload)

    except Exception as exc:
        error_response(
            "Failed to query Google Analytics Data API",
            details={
                "message": str(exc),
                "traceback": traceback.format_exc().splitlines(),
            },
        )


if __name__ == "__main__":
    main()
