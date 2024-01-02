document.addEventListener('DOMContentLoaded', () => {
    populateGeneCartDropdown();
});

const populateGeneCartDropdown = async () => {
    let numEntries = 10;

    try {
        const data = await apiCallsMixin.fetchGeneCarts();
        console.log(data);
        const template = document.querySelector('#tmpl-gene-list-item');

        for (const entry of data.domain_carts) {
            const row = template.content.cloneNode(true);
            row.querySelector('.gene-list-item-label').textContent = entry.label;
            //document.querySelector('#dropdown-content-gene-lists').appendChild(row);
            // reduce number of entries by 1
            if (numEntries > 0) {
                document.querySelector('#dropdown-content-gene-lists').appendChild(row);
                numEntries--;
            }
        }

    } catch (error) {
        console.error(error);
    }
}

const populateUserHistoryTable = async () => {
    const numEntries = 5;

    // Load the spinner template
    const spinnerTemplate = document.querySelector('#user-history-loading');
    document.querySelector('#user-history-table-tbody').innerHTML = '';
    document.querySelector('#user-history-table-tbody').appendChild(spinnerTemplate.content.cloneNode(true));

    try {
        const data = await apiCallsMixin.fetchUserHistoryEntries(numEntries);
        const template = document.querySelector('#user-history-row');
        document.querySelector('#user-history-table-tbody').innerHTML = '';

        if (data.length === 0) {
            const noHistoryTemplate = document.querySelector('#user-history-no-entries');
            document.querySelector('#user-history-table-tbody').appendChild(noHistoryTemplate.content.cloneNode(true));
        } else {
            for (const entry of data) {
                const row = template.content.cloneNode(true);

                let formatted_category = entry.entry_category.replaceAll('_', ' ');
                formatted_category = formatted_category.charAt(0).toUpperCase() + formatted_category.slice(1);

                row.querySelector('.category').textContent = formatted_category;
                row.querySelector('.action-label').textContent = entry.label;
                row.querySelector('.date').textContent = entry.entry_date;
                row.querySelector('.url').setAttribute('href', entry.url);

                document.querySelector('#user-history-table-tbody').appendChild(row);
            }
        }
    } catch (error) {
        console.error(error);
    }
}

const handlePageSpecificLoginUIUpdates = async (event) => {
    if (CURRENT_USER.session_id) {
        populateUserHistoryTable();
    }
}
