document.addEventListener('DOMContentLoaded', () => {
    // Nothing here yet
});

const populateUserHistoryTable = async () => {
    const numEntries = 5;
    try {
        const data = await apiCallsMixin.fetchUserHistoryEntries(numEntries);
        const template = document.querySelector('#user-history-row');

        for (const entry of data) {
            const row = template.content.cloneNode(true);

            let formatted_category = entry.entry_category.replaceAll('_', ' ');
            formatted_category = formatted_category.charAt(0).toUpperCase() + formatted_category.slice(1);

            row.querySelector('.category').textContent = formatted_category;
            row.querySelector('.action-label').textContent = entry.label;
            row.querySelector('.date').textContent = entry.entry_date;
            row.querySelector('.url').setAttribute('href', entry.url);

            document.querySelector('#user-history-table-tbody').appendChild(row);
        };
    } catch (error) {
        console.error(error);
    }
}

const handlePageSpecificLoginUIUpdates = async (event) => {
    if (CURRENT_USER) {
        populateUserHistoryTable();
    }
}
