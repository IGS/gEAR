window.onload=function() {
    // check if the user is already logged in
    check_for_login();
    session_id = Cookies.get('gear_session_id');
}
