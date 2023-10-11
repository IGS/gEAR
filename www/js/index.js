document.addEventListener('DOMContentLoaded', () => {
    if (CURRENT_USER) {
        populateUserHistoryTable();
    }
});

const populateUserHistoryTable = async (el, eventType) => {
    const payload = { 'session_id': CURRENT_USER.session_id,
                      'entries': 5
                    }
    const {data} = await axios.post("/cgi/get_user_history_entries.cgi", convertToFormData(payload));
    console.log(data);
}

const handlePageSpecificLoginUIUpdates = async (event) => {
    // nothing here yet
}
