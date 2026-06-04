from google.analytics.data_v1beta import BetaAnalyticsDataClient
from google.analytics.data_v1beta.types import DateRange, Dimension, Metric, RunReportRequest


def run_report(property_id: str):
    client = BetaAnalyticsDataClient()

    request = RunReportRequest(
        property=f"properties/{property_id}",
        dimensions=[
            Dimension(name="date"),
        ],
        metrics=[
            Metric(name="activeUsers"),
            Metric(name="sessions"),
            Metric(name="screenPageViews"),
        ],
        date_ranges=[
            DateRange(start_date="30daysAgo", end_date="today"),
        ],
    )

    response = client.run_report(request)

    rows = []
    for row in response.rows:
        rows.append({
            "date": row.dimension_values[0].value,
            "activeUsers": row.metric_values[0].value,
            "sessions": row.metric_values[1].value,
            "screenPageViews": row.metric_values[2].value,
        })

    return rows


if __name__ == "__main__":
    PROPERTY_ID = "312199300"
    data = run_report(PROPERTY_ID)
    for row in data:
        print(row)