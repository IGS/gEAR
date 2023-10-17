document.addEventListener('DOMContentLoaded', () => {
    // Nothing here yet
});

const populateUserHistoryTable = async (el, eventType) => {
    const payload = { 'session_id': CURRENT_USER.session_id,
                      'entries': 5
                    }
    const {data} = await axios.post("/cgi/get_user_history_entries.cgi", convertToFormData(payload));
    const template = document.querySelector('#user-history-row');

    data.forEach(function (entry, idx) {
        const row = template.content.cloneNode(true);

        let formatted_category = entry.entry_category.replaceAll('_', ' ');
        formatted_category = formatted_category.charAt(0).toUpperCase() + formatted_category.slice(1);
        
        row.querySelector('.category').textContent = formatted_category;
        row.querySelector('.action-label').textContent = entry.label;
        row.querySelector('.date').textContent = entry.entry_date;
        document.querySelector('#user-history-table-tbody').appendChild(row);
    });
}

const handlePageSpecificLoginUIUpdates = async (event) => {
    if (CURRENT_USER) {
        populateUserHistoryTable();
    }
}
