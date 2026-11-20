/*
 * gEAR analytics dashboard.
 *
 * Data sources:
 *   - Site traffic:     /cgi/stats.cgi (Google Analytics 4)
 *   - Platform content: /cgi/get_dataset_analytics.cgi (database)
 */

/* --------- */
/* Constants */
/* --------- */
const ANALYTICS_URL = '/cgi/stats.cgi';
const DATASET_STATS_URL = '/cgi/get_dataset_analytics.cgi';
const ANALYTICS_DAYS = 30;
const RANGE_PRESETS = [7, 30, 90, 365];
const TOP_N = 10;
const BYTES_PER_GB = 1024 ** 3;
const MS_PER_DAY = 86400000;
const LIVE_BLUE = '#2f7ee8';

const COLORS = {
    primary: '#362E5C',
    secondary: '#544A8E',
    teal: '#1b9aaa',
    orange: '#e08a1e',
    grid: '#ebe9f4',
    text: '#5c5877'
};
const PALETTE = ['#544A8E', '#1b9aaa', '#e08a1e', '#c8556e', '#7bb369', '#8d99ae', '#362E5C', '#b08ad6'];

const charts = {};

const state = {
    // { mode: 'preset', days } or { mode: 'custom', start: 'YYYY-MM-DD', end: 'YYYY-MM-DD' }
    range: { mode: 'preset', days: ANALYTICS_DAYS },
    results: { analytics: null, datasets: null },
    pending: 0
};

const tables = new Map();

/* ------- */
/* Helpers */
/* ------- */
const escapeHtml = (value) => String(value ?? '')
    .replace(/&/g, '&amp;')
    .replace(/</g, '&lt;')
    .replace(/>/g, '&gt;')
    .replace(/"/g, '&quot;')
    .replace(/'/g, '&#39;');

const formatNumber = (num) => Number(num ?? 0).toLocaleString('en-US');

// stats.cgi returns dates as YYYYMMDD
const parseAnalyticsDate = (dateStr) => new Date(
    Number(dateStr.substring(0, 4)),
    Number(dateStr.substring(4, 6)) - 1,
    Number(dateStr.substring(6, 8))
);

const formatAnalyticsDate = (dateStr) => parseAnalyticsDate(dateStr)
    .toLocaleDateString('en-US', { month: 'short', day: 'numeric' });

const formatFileSize = (sizeInGB) => {
    const units = ['B', 'KB', 'MB', 'GB', 'TB', 'PB'];
    let size = Number(sizeInGB ?? 0) * BYTES_PER_GB;
    let unitIndex = 0;
    while (size >= 1024 && unitIndex < units.length - 1) {
        size /= 1024;
        unitIndex++;
    }
    return `${formatNumber(+size.toFixed(2))} ${units[unitIndex]}`;
};

const capitalizeFirstLetter = (str) => (str ? str.charAt(0).toUpperCase() + str.slice(1) : '');

const getById = (id) => document.getElementById(id);

const setHtml = (id, html) => {
    const el = getById(id);
    if (el) el.innerHTML = html;
};

const renderEmpty = (message) => `<p class="empty-state">${escapeHtml(message)}</p>`;

const renderChart = (key, canvasId, config) => {
    const canvas = getById(canvasId);
    if (!canvas || typeof Chart === 'undefined') return;
    charts[key]?.destroy();
    charts[key] = new Chart(canvas.getContext('2d'), config);
};

const doughnutConfig = (labels, data) => ({
    type: 'doughnut',
    data: {
        labels,
        datasets: [{
            data,
            backgroundColor: labels.map((_, i) => PALETTE[i % PALETTE.length]),
            borderColor: '#fff',
            borderWidth: 2
        }]
    },
    options: {
        responsive: true,
        maintainAspectRatio: false,
        cutout: '60%',
        plugins: {
            legend: { position: 'bottom', labels: { boxWidth: 12, color: COLORS.text } }
        }
    }
});

/* ------------ */
/* Data loading */
/* ------------ */
const fetchJson = async (url) => {
    const response = await fetch(url);
    // CGI error payloads are JSON even on non-200 statuses
    const data = await response.json().catch(() => null);
    if (!data) throw new Error(`Unexpected response from ${url} (HTTP ${response.status})`);
    return data;
};

/* ---- Date range helpers ---- */
const pad2 = (n) => String(n).padStart(2, '0');
const toIsoDate = (date) => `${date.getFullYear()}-${pad2(date.getMonth() + 1)}-${pad2(date.getDate())}`;
const toAnalyticsDate = (date) => toIsoDate(date).replace(/-/g, '');
const parseIsoDate = (iso) => {
    const [year, month, day] = iso.split('-').map(Number);
    return new Date(year, month - 1, day);
};
const addDays = (date, days) => new Date(date.getFullYear(), date.getMonth(), date.getDate() + days);

const formatLongDate = (date) => date.toLocaleDateString('en-US', { month: 'short', day: 'numeric', year: 'numeric' });

const describeRange = (range) => (range.mode === 'custom'
    ? `${formatLongDate(parseIsoDate(range.start))} – ${formatLongDate(parseIsoDate(range.end))}`
    : `Last ${range.days} days`);

const buildAnalyticsUrl = (range) => {
    const params = new URLSearchParams({ top_n: TOP_N });
    if (range.mode === 'custom') {
        params.set('start_date', range.start);
        params.set('end_date', range.end);
    } else {
        params.set('days', range.days);
    }
    return `${ANALYTICS_URL}?${params}`;
};

const loadAnalytics = async (range) => {
    const data = await fetchJson(buildAnalyticsUrl(range));
    if (!data.ok) throw new Error(data.error || 'Google Analytics request failed');
    return { data, isTest: false };
};

const loadDatasetStats = async () => {
    const data = await fetchJson(DATASET_STATS_URL);
    if (data.error) throw new Error(data.error);
    return { data, isTest: false };
};

const cleanupData = (data) => {
    if (!data) return;
    if (Array.isArray(data.topPages)) {
        const cleaned = { };
        for (const page of data.topPages) {
            let path = page.pagePath || '/';
            switch (page.pagePath) {
                case '/':
                    path = '/index.html';
                    break;
            }

            if (path in cleaned) {
                // Merge metrics
                const target = cleaned[path];
                target.screenPageViews = (target.screenPageViews || 0) + (page.screenPageViews || 0);
                target.activeUsers = (target.activeUsers || 0) + (page.activeUsers || 0);
                target.sessions = (target.sessions || 0) + (page.sessions || 0);
                target.allPagePaths += ` · ${page.pagePath}`;
            } else {
                cleaned[path] = page;
                page.allPagePaths = page.pagePath;
                page.pagePath = path;
            }
        }

        data.topPages = Object.values(cleaned);
    }
};

/* ----------------- */
/* Status & errors */
/* ----------------- */
const clearError = (key) => getById(`error-${key}`)?.remove();

const showError = (key, title, error) => {
    console.error(title, error);
    clearError(key);
    getById('errorRegion')?.insertAdjacentHTML('beforeend', `
        <div id="error-${key}" class="notification is-danger is-light">
            <strong>${escapeHtml(title)}</strong> ${escapeHtml(error?.message ?? error)}
        </div>
    `);
};

const setLoading = (delta) => {
    state.pending += delta;
    const isLoading = state.pending > 0;
    getById('loadingBar')?.classList.toggle('is-hidden', !isLoading);
    const button = getById('refreshBtn');
    if (button) button.disabled = isLoading;
};

const updateSourceTag = () => {
    const tag = getById('dataSourceTag');
    if (!tag) return;
    const loaded = Object.values(state.results).filter(Boolean);
    const testCount = loaded.filter((r) => r.isTest).length;

    tag.classList.remove('is-warning', 'is-light');
    if (!loaded.length) {
        tag.classList.add('is-light');
        tag.textContent = 'No data';
    } else if (testCount === loaded.length) {
        tag.classList.add('is-warning');
        tag.textContent = 'Test data';
    } else if (testCount > 0) {
        tag.classList.add('is-warning');
        tag.textContent = 'Mixed test and live data';
    } else {
        tag.classList.add('is-light');
        tag.textContent = 'Live data';
    }
};

/* --------- */
/* Renderers */
/* --------- */

/* ---- Site traffic renderers ---- */
const renderKpi = ({ icon, value, label }) => `
    <div class="kpi">
        <span class="kpi-icon" aria-hidden="true"><i class="mdi ${icon} mdi-24px"></i></span>
        <div>
            <div class="kpi-value">${escapeHtml(value)}</div>
            <div class="kpi-label">${escapeHtml(label)}</div>
        </div>
    </div>
`;

const renderSummaryCards = (summary) => {
    const engagementRate = summary.sessions ? Math.round((summary.engagedSessions / summary.sessions) * 100) : 0;
    const cards = [
        { icon: 'mdi-account-group-outline', value: formatNumber(summary.activeUsers), label: 'Active users' },
        { icon: 'mdi-clock-outline', value: formatNumber(summary.sessions), label: 'Sessions' },
        { icon: 'mdi-file-eye-outline', value: formatNumber(summary.screenPageViews), label: 'Page views' },
        { icon: 'mdi-account-heart-outline', value: formatNumber(summary.engagedSessions), label: 'Engaged sessions' },
        { icon: 'mdi-percent-outline', value: `${engagementRate}%`, label: 'Engagement rate' },
        { icon: 'mdi-cursor-default-click-outline', value: formatNumber(summary.eventCount), label: 'Total events' }
    ];
    setHtml('summaryCards', cards.map(renderKpi).join(''));
};

const renderTimeseriesChart = (timeseries) => {
    const sorted = [...timeseries].sort((a, b) => a.date.localeCompare(b.date));
    const line = (label, key, color, yAxisID) => ({
        label,
        data: sorted.map((d) => d[key]),
        borderColor: color,
        backgroundColor: color,
        borderWidth: 2,
        pointRadius: 0,
        pointHoverRadius: 4,
        tension: 0.3,
        yAxisID
    });
    const axisTitle = (text) => ({ display: true, text, color: COLORS.text });

    renderChart('timeseries', 'timeseriesChart', {
        type: 'line',
        data: {
            labels: sorted.map((d) => formatAnalyticsDate(d.date)),
            datasets: [
                line('Active users', 'activeUsers', COLORS.primary, 'y'),
                line('Sessions', 'sessions', COLORS.teal, 'y'),
                line('Page views', 'screenPageViews', COLORS.orange, 'y1')
            ]
        },
        options: {
            responsive: true,
            maintainAspectRatio: false,
            interaction: { mode: 'index', intersect: false },
            plugins: { legend: { position: 'top', align: 'end', labels: { boxWidth: 12, color: COLORS.text } } },
            scales: {
                x: { grid: { display: false }, ticks: { maxTicksLimit: 10, color: COLORS.text } },
                y: {
                    beginAtZero: true,
                    position: 'left',
                    title: axisTitle('Users / sessions'),
                    grid: { color: COLORS.grid },
                    ticks: { color: COLORS.text }
                },
                y1: {
                    beginAtZero: true,
                    position: 'right',
                    title: axisTitle('Page views'),
                    grid: { drawOnChartArea: false },
                    ticks: { color: COLORS.text }
                }
            }
        }
    });
};

const renderTable = (containerId, caption, columns, rows) => {
    if (!tables.has(containerId)) {
        tables.set(containerId, { sortIndex: columns.findIndex((col) => col.kind === 'number'), direction: 'desc' });
    }
    Object.assign(tables.get(containerId), { caption, columns, rows: rows ?? [] });
    drawTable(containerId);
};

const drawTable = (containerId) => {
    const { caption, columns, rows, sortIndex, direction } = tables.get(containerId);

    if (!rows.length) {
        setHtml(containerId, renderEmpty('No data for this period.'));
        return;
    }

    const sortValue = columns[sortIndex].value;
    const sign = direction === 'asc' ? 1 : -1;
    const sorted = [...rows].sort((a, b) => sign * (sortValue(a) - sortValue(b)));
    const barMax = Math.max(...rows.map(sortValue), 1);

    const renderHeader = (col, index) => {
        if (col.kind !== 'number') return `<th scope="col">${escapeHtml(col.header)}</th>`;
        const isSorted = index === sortIndex;
        const icon = !isSorted ? 'mdi-unfold-more-horizontal' : (direction === 'asc' ? 'mdi-arrow-up' : 'mdi-arrow-down');
        const ariaSort = !isSorted ? 'none' : (direction === 'asc' ? 'ascending' : 'descending');
        return `<th scope="col" class="num" aria-sort="${ariaSort}">
            <button type="button" class="sort-btn${isSorted ? ' is-active' : ''}" data-col="${index}">
                <span>${escapeHtml(col.header)}</span>
                <span class="icon is-small" aria-hidden="true"><i class="mdi ${icon}"></i></span>
            </button>
        </th>`;
    };

    const renderCell = (col, row, index) => {
        const value = col.value(row);
        const title = col.title?.(row) ?? value;
        if (col.kind === 'number') {
            if (index !== sortIndex) return `<td class="num">${formatNumber(value)}</td>`;
            return `<td class="num"><div class="bar-cell">
                <span>${formatNumber(value)}</span>
                <span class="bar-track" aria-hidden="true"><span class="bar-fill" style="width:${Math.round((value / barMax) * 100)}%"></span></span>
            </div></td>`;
        }
        const text = col.kind === 'code' ? `<code>${escapeHtml(value)}</code>` : escapeHtml(value);
        return `<td class="label-cell" title="${escapeHtml(title)}">${text}</td>`;
    };

    setHtml(containerId, `
        <table class="dash-table">
            <caption class="is-sr-only">${escapeHtml(caption)}</caption>
            <thead>
                <tr>${columns.map(renderHeader).join('')}</tr>
            </thead>
            <tbody>
                ${sorted.map((row) => `<tr>${columns.map((col, i) => renderCell(col, row, i)).join('')}</tr>`).join('')}
            </tbody>
        </table>
    `);
};

const handleSortClick = (event) => {
    const button = event.target.closest('.sort-btn');
    const container = button?.closest('[data-sortable-table]');
    const table = tables.get(container?.id);
    if (!table) return;

    const index = Number(button.dataset.col);
    if (index === table.sortIndex) {
        table.direction = table.direction === 'desc' ? 'asc' : 'desc';
    } else {
        table.sortIndex = index;
        table.direction = 'desc';
    }
    drawTable(container.id);
    // drawTable replaces the button, so restore keyboard focus
    container.querySelector(`.sort-btn[data-col="${index}"]`)?.focus();
};

const renderRealtime = (realtime) => {
    if (!realtime) {
        setHtml('realtimeContainer', renderEmpty('Real-time data is unavailable.'));
        return;
    }

    const total = realtime.activeUsersNow ?? 0;
    const countries = [...(realtime.byCountry ?? [])].sort((a, b) => b.activeUsers - a.activeUsers);
    const max = Math.max(...countries.map((c) => c.activeUsers), 1);

    const list = countries.map((c) => `
        <li>
            <span class="live-country" title="${escapeHtml(c.country)}">${escapeHtml(c.country)}</span>
            <span class="live-bar" aria-hidden="true"><span style="width:${Math.round((c.activeUsers / max) * 100)}%"></span></span>
            <strong>${formatNumber(c.activeUsers)}</strong>
        </li>`).join('');

    setHtml('realtimeContainer', `
        <div class="live-hero">
            <div class="live-number">${formatNumber(total)}</div>
            <div class="live-caption">active user${total === 1 ? '' : 's'} right now</div>
        </div>
        ${countries.length ? `<h4 class="live-subtitle">By country</h4><ul class="live-list">${list}</ul>` : ''}
    `);
};

const renderAnalytics = (data) => {
    const rangeLabel = getById('trafficRangeLabel');
    if (rangeLabel) rangeLabel.textContent = `${describeRange(state.range)}, from Google Analytics`;

    cleanupData(data); // some data cleanup

    renderSummaryCards(data.summary);
    renderTimeseriesChart(data.timeseries ?? []);
    renderRealtime(data.realtime);

    renderTable('topPagesContainer', 'Top pages', [
        { header: 'Page', value: (r) => r.pagePath, title: (r) => r.allPagePaths, kind: 'code' },
        { header: 'Page views', value: (r) => r.screenPageViews, kind: 'number' },
        { header: 'Users', value: (r) => r.activeUsers, kind: 'number' },
        { header: 'Sessions', value: (r) => r.sessions, kind: 'number' }
    ], data.topPages);

    renderTable('trafficSourcesContainer', 'Traffic sources', [
        { header: 'Source / medium', value: (r) => r.sessionSourceMedium },
        { header: 'Sessions', value: (r) => r.sessions, kind: 'number' },
        { header: 'Users', value: (r) => r.activeUsers, kind: 'number' }
    ], data.trafficSources);

    renderTable('countriesContainer', 'Countries', [
        { header: 'Country', value: (r) => r.country },
        { header: 'Users', value: (r) => r.activeUsers, kind: 'number' },
        { header: 'Sessions', value: (r) => r.sessions, kind: 'number' },
        { header: 'Page views', value: (r) => r.screenPageViews, kind: 'number' }
    ], data.countries);

    renderTable('topEventsContainer', 'Top events', [
        { header: 'Event', value: (r) => r.eventName, kind: 'code' },
        { header: 'Count', value: (r) => r.eventCount, kind: 'number' },
        { header: 'Users', value: (r) => r.activeUsers, kind: 'number' }
    ], data.topEvents);
};

/* ---- Platform content renderers ---- */
const renderPlatformCards = (stats) => {
    const overview = stats.overview ?? {};
    const statistics = stats.statistics ?? {};
    const cards = [
        { icon: 'mdi-database-outline', value: formatNumber(overview.numberOfDatasets), label: 'Total datasets' },
        { icon: 'mdi-web', value: formatNumber(overview.publicDatasets), label: 'Public datasets' },
        { icon: 'mdi-lock-outline', value: formatNumber(overview.privateDatasets), label: 'Private datasets' },
        { icon: 'mdi-file-send-outline', value: formatNumber(statistics.sharedDatasets), label: 'Shared datasets' },
        { icon: 'mdi-content-save-outline', value: formatFileSize(overview.totalSize), label: 'Total storage' },
        { icon: 'mdi-file-chart-outline', value: formatFileSize(overview.averageDatasetSize), label: 'Avg. dataset size' },
        { icon: 'mdi-account-multiple-outline', value: formatNumber(statistics.registeredUsers), label: 'Registered users' },
        { icon: 'mdi-share-variant-outline', value: formatNumber(statistics.totalShares), label: 'Total shares' },
        { icon: 'mdi-dna', value: formatNumber(statistics.totalOrganisms), label: 'Organisms' },
        { icon: 'mdi-format-list-bulleted', value: formatNumber(statistics.geneLists), label: 'Gene lists' },
        { icon: 'mdi-view-dashboard-outline', value: formatNumber(statistics.layouts), label: 'Layouts' }
    ];
    setHtml('platformCards', cards.map(renderKpi).join(''));
};

// Cumulative dataset count per month, with empty months filled in
const buildCumulativeMonthly = (datesAdded) => {
    const perMonth = {};
    for (const [date, count] of Object.entries(datesAdded)) {
        const month = date.substring(0, 7);
        perMonth[month] = (perMonth[month] ?? 0) + count;
    }

    const months = Object.keys(perMonth).sort();
    if (!months.length) return { labels: [], data: [] };

    const [lastYear, lastMonth] = months[months.length - 1].split('-').map(Number);
    let [year, month] = months[0].split('-').map(Number);
    const labels = [];
    const data = [];
    let running = 0;
    while (year < lastYear || (year === lastYear && month <= lastMonth)) {
        const key = `${year}-${String(month).padStart(2, '0')}`;
        running += perMonth[key] ?? 0;
        labels.push(new Date(year, month - 1, 1).toLocaleDateString('en-US', { month: 'short', year: 'numeric' }));
        data.push(running);
        month++;
        if (month > 12) {
            month = 1;
            year++;
        }
    }
    return { labels, data };
};

const renderDatasetCharts = (stats) => {
    const { labels, data } = buildCumulativeMonthly(stats.statistics?.datesAdded ?? {});
    renderChart('datasetsDates', 'datasetsDatesChart', {
        type: 'line',
        data: {
            labels,
            datasets: [{
                label: 'Total datasets',
                data,
                borderColor: COLORS.secondary,
                backgroundColor: 'rgba(84, 74, 142, 0.15)',
                borderWidth: 2,
                fill: true,
                tension: 0.3,
                pointRadius: 0,
                pointHoverRadius: 4
            }]
        },
        options: {
            responsive: true,
            maintainAspectRatio: false,
            interaction: { mode: 'index', intersect: false },
            plugins: { legend: { display: false } },
            scales: {
                x: { grid: { display: false }, ticks: { maxTicksLimit: 8, color: COLORS.text } },
                y: { beginAtZero: true, grid: { color: COLORS.grid }, ticks: { color: COLORS.text } }
            }
        }
    });

    // byType keys are already friendly names (see get_dataset_analytics.cgi)
    const byType = stats.byType ?? {};
    const typeLabels = Object.keys(byType);
    renderChart('datasetsType', 'datasetsTypeChart', doughnutConfig(typeLabels, typeLabels.map((k) => byType[k].count)));

    const status = stats.statusBreakdown ?? {};
    const statusKeys = Object.keys(status);
    renderChart('datasetsStatus', 'datasetsStatusChart', doughnutConfig(statusKeys.map(capitalizeFirstLetter), statusKeys.map((k) => status[k])));
};

const renderDatasetStats = (stats) => {
    renderPlatformCards(stats);
    renderDatasetCharts(stats);
};

/* ---- Date range controls ---- */
let customRangeOpen = false;

const syncRangeControls = () => {
    const { range } = state;
    for (const button of document.querySelectorAll('[data-range-days]')) {
        const isActive = range.mode === 'preset' && Number(button.dataset.rangeDays) === range.days;
        button.classList.toggle('is-active', isActive);
        button.setAttribute('aria-pressed', String(isActive));
    }

    const customButton = getById('rangeCustomBtn');
    const customActive = range.mode === 'custom' || customRangeOpen;
    customButton?.classList.toggle('is-active', customActive);
    customButton?.setAttribute('aria-pressed', String(customActive));
    getById('customRangeForm')?.classList.toggle('is-hidden', !customActive);

    const start = getById('rangeStart');
    const end = getById('rangeEnd');
    const today = toIsoDate(new Date());
    if (start && end) {
        start.max = today;
        end.max = today;
        if (range.mode === 'custom') {
            start.value = range.start;
            end.value = range.end;
        } else if (!start.value && !end.value) {
            start.value = toIsoDate(addDays(new Date(), -(range.days - 1)));
            end.value = today;
        }
    }
};

const setRangeError = (message) => {
    const el = getById('rangeError');
    if (el) el.textContent = message;
};

const selectPresetRange = (days) => {
    customRangeOpen = false;
    setRangeError('');
    state.range = { mode: 'preset', days };
    syncRangeControls();
    loadAnalyticsSection();
};

const applyCustomRange = (event) => {
    event.preventDefault();
    const start = getById('rangeStart').value;
    const end = getById('rangeEnd').value;

    if (!start || !end) {
        setRangeError('Choose both a start and an end date.');
    } else if (start > end) {
        setRangeError('The start date must be on or before the end date.');
    } else if (end > toIsoDate(new Date())) {
        setRangeError('The end date cannot be in the future.');
    } else {
        setRangeError('');
        customRangeOpen = false;
        state.range = { mode: 'custom', start, end };
        syncRangeControls();
        loadAnalyticsSection();
    }
};

/* -------------- */
/* Initialization */
/* -------------- */
const markUpdated = () => {
    const lastUpdated = getById('lastUpdated');
    if (lastUpdated) lastUpdated.textContent = `Updated ${new Date().toLocaleTimeString('en-US')}`;
};

// Guards against an older, slower response overwriting a newer one
let latestAnalyticsRequest = 0;
const loadAnalyticsSection = async () => {
    const requestId = ++latestAnalyticsRequest;
    setLoading(1);
    try {
        const result = await loadAnalytics(state.range);
        if (requestId !== latestAnalyticsRequest) return;
        state.results.analytics = result;
        clearError('analytics');
        renderAnalytics(result.data);
    } catch (error) {
        if (requestId !== latestAnalyticsRequest) return;
        state.results.analytics = null;
        showError('analytics', 'Could not load site traffic.', error);
    } finally {
        setLoading(-1);
        updateSourceTag();
        markUpdated();
    }
};

const loadDatasetSection = async () => {
    setLoading(1);
    try {
        const result = await loadDatasetStats();
        state.results.datasets = result;
        clearError('datasets');
        renderDatasetStats(result.data);
    } catch (error) {
        state.results.datasets = null;
        showError('datasets', 'Could not load platform content.', error);
    } finally {
        setLoading(-1);
        updateSourceTag();
        markUpdated();
    }
};

const loadDashboard = () => {
    if (typeof Chart === 'undefined') {
        showError('charts', 'Charts unavailable.', 'The Chart.js library could not be loaded.');
    }
    return Promise.all([loadAnalyticsSection(), loadDatasetSection()]);
};

/* --------------- */
/* Event listeners */
/* --------------- */
getById('refreshBtn')?.addEventListener('click', loadDashboard);

for (const button of document.querySelectorAll('[data-range-days]')) {
    button.addEventListener('click', () => selectPresetRange(Number(button.dataset.rangeDays)));
}

getById('rangeCustomBtn')?.addEventListener('click', () => {
    customRangeOpen = true;
    syncRangeControls();
    getById('rangeStart')?.focus();
});

getById('customRangeForm')?.addEventListener('submit', applyCustomRange);

for (const container of document.querySelectorAll('[data-sortable-table]')) {
    container.addEventListener('click', handleSortClick);
}

syncRangeControls();
loadDashboard();
