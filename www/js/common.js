let CURRENT_USER = null;
let SITE_PREFS = null;

let session_id = null;

// extends JQuery to provide an .exists() function for any selector
$.fn.exists = function () {
    return this.length !== 0;
}

const sleep = (milliseconds) => {
  return new Promise(resolve => setTimeout(resolve, milliseconds))
}

// Exclude SVG path elements. We generate their tooltips in main.js:select_search_result()
const elementsWithTipsList = '[title]:not("path")';

// Do any includes first
$(document).ready(() => {
    $('#navigation_bar').load('./include/navigation_bar.html', () => {
        // Load popover info
        load_forgot_password();

        // Now that the navigation bar is loaded, we can check if the user is logged in
        // and by extension, handle the login UI updates
        check_for_login();

        if (! show_video_link()) {
            $(".js-video-link").hide();
        }
    });

    $.ajax({
        url: "/site_domain_prefs.json",
        dataType: 'json',
        //async: false,
     }).done((data) => {
        SITE_PREFS = data;
        loadCSS(`./css/by_domain/${SITE_PREFS['domain_label']}/theme_colors.${(new Date()).getTime()}.css`);
        $('#funding').load(`./include/by_domain/${SITE_PREFS['domain_label']}/funding.html`);
        $('#footer').load(`./include/by_domain/${SITE_PREFS['domain_label']}/footer.html`);
        $('#site_label_c').load(`./include/by_domain/${SITE_PREFS['domain_label']}/site_label_bar.html`);
        $('a#main_logo').css('background-image', `url("../img/by_domain/${SITE_PREFS['domain_label']}/logo_standard.png")`);

        const title_page_root = `./include/by_domain/${SITE_PREFS['domain_label']}/page_title_root.html`;
        $.get(title_page_root, (data) => {
            $('title').html(`${data} - ${$('title').html()}`);
        });

        // now check and make sure these steps happened, and try again if not.  WHY???
        // Every time this happens the steps after funding loader fail.  Retry.
        sleep(500).then(() => {
            // did the background image load?
            bg = $('a#main_logo').css('background-image');
            if (bg == 'none') {
                $('#footer').load(`./include/by_domain/${SITE_PREFS['domain_label']}/footer.html`);
                $('#site_label_c').load(`./include/by_domain/${SITE_PREFS['domain_label']}/site_label_bar.html`);
                $('a#main_logo').css('background-image', `url("../img/by_domain/${SITE_PREFS['domain_label']}/logo_standard.png")`);
            }
        });

        // populate any site-specific labels, usually spans
        $('.domain_short_display_label').text(SITE_PREFS['domain_short_display_label']);

        let page_name = location.pathname;
        if (page_name == "/") {
            page_name = 'index.html';
        }

        const head = document.getElementsByTagName('head')[0];
        const body = document.getElementsByTagName('body')[0];

        // load analytics
        const ga_script = document.createElement('script');
        //ga_script.src = "https://www.google-analytics.com/analytics.js";
        ga_script.src = "https://www.googletagmanager.com/gtag/js?id=" + SITE_PREFS['google_analytics_4_measurement_id'];
        ga_script.async = ""
        head.append(ga_script)

        window.dataLayer = window.dataLayer || [];
        function gtag(){dataLayer.push(arguments);}
        gtag('js', new Date());
        gtag('config', SITE_PREFS['google_analytics_4_measurement_id']);

        // Load plugins, if any
        for (const [plugin_name, plugin_page_names] of Object.entries(SITE_PREFS['enabled_plugins'])) {
            if (plugin_page_names.includes(page_name)) {
                var plugin_import_basename = page_name.replace(/\.[^/.]+$/, "");

                // create a hidden element at the end of the document and put it there
                var plugin_import_html_url = "./plugins/" + plugin_name + "/" + page_name;
                var plugin_html_element = document.createElement('div');
                plugin_html_element.id = plugin_name + "_html_c";
                plugin_html_element.style = "display: none;"
                body.append(plugin_html_element)
                $.get(plugin_import_html_url, (data) => {
                    $(`#${plugin_name}_html_c`).html(data);
                });

                var plugin_import_css_url = "./plugins/" + plugin_name + "/" + plugin_import_basename + ".css";
                var style = document.createElement('link');
                style.href = plugin_import_css_url;
                style.type = 'text/css';
                style.rel = 'stylesheet';
                head.append(style);

                var plugin_import_js_url = "./plugins/" + plugin_name + "/" + plugin_import_basename + ".js";
                var script = document.createElement('script');
                script.src = plugin_import_js_url;
                script.type = 'text/javascript';
                head.append(script);
            }
        }
    }).fail((jqXHR, textStatus, errorThrown) => {
        console.error(`Error loading site_domain_prefs.json: ${textStatus}`);
    });

    check_browser();

});

function check_browser() {
    /*
    Checks if user is using Chrome browser. If not Chrome, alerts user.
    Sets a cookie for so user will not get alerted after initial alert. Expires 1 day.
    */

    let isChrome = Cookies.get('gear_browser_ischrome');
    if (typeof isChrome === 'undefined') {
        // Check if the browser is Chrome
        isChrome = navigator.userAgent.toLowerCase().includes('chrome');
        if (!isChrome) {
            //Alert user the browser is not supported.
            $('.alert-container').html('<div class="alert alert-danger alert-dismissible" role="alert">' +
                  '<button type="button" class="close close-alert" data-dismiss="alert" aria-label="Close"><span aria-hidden="true">&times;</span></button>' +
                  '<p class="alert-message">' +
                  'The gEAR is not optimized for this browser, so its appearance and functions may not function correctly.</p>' +
                  '<p class="alert-message">' +
                  'Please use the Chrome browser. If you do not have Chrome, please download it <a href="https://www.google.com/chrome/" target="blank" title="Link to download Chrome">here</a>.' +
                  '</p>' +
                  '</div>').show();
        }

        //Set cookie so user does not repeatedly get Unsupported browser message
        Cookies.set('gear_browser_ischrome', isChrome, { expires: 1 });
    } else if (!isChrome) {
        console.warn("Unsupported browser detected");
    }
}; //end check_browser()

function check_for_login() {
    // success:  returns session_id
    // failure:  returns null
    session_id = Cookies.get('gear_session_id');
    CURRENT_USER = new User();

    if (typeof session_id !== 'undefined') {
        // we found a cookie, so see if it's a valid session
        $.ajax({
            url : './cgi/get_session_info.cgi',
            type: "POST",
            async: false,
            data : { 'session_id': session_id },
            dataType:"json"
        }).done((data) => {
            CURRENT_USER = new User({session_id, ...data});
            CURRENT_USER.set_default_profile();
            $('span.user_logged_in').text(CURRENT_USER.user_name);
            handle_login_ui_updates();
        }).fail((jqXHR, textStatus, errorThrown) => {
            console.error(`Error getting session info: ${textStatus}`);
        });
    } else {
        handle_login_ui_updates();
    }

    // Now that session_id has been obtained, we can trigger events that depend on it
    $(document).trigger("build_jstrees");

    return session_id;
}

function common_datetime() {
    const today = new Date();
    const date = `${today.getFullYear()}-${today.getMonth()+1}-${today.getDate()}`;
    const time = `${today.getHours()}:${today.getMinutes()}:${today.getSeconds()}`;
    return `${date} ${time}`;
}

function copyToClipboard(text) {
    // https://stackoverflow.com/a/59594066
    if (window.clipboardData?.setData) {
        // IE specific code path to prevent textarea being shown while dialog is visible.
        return clipboardData.setData("Text", text);

    } else if (document.queryCommandSupported?.("copy")) {
        const textarea = document.createElement("textarea");
        textarea.textContent = text;
        textarea.style.position = "fixed";  // Prevent scrolling to bottom of page in MS Edge.
        document.body.appendChild(textarea);
        textarea.select();
        try {
            return document.execCommand("copy");  // Security exception may be thrown by some browsers.
        } catch (ex) {
            console.warn("Copy to clipboard failed.", ex);
            return false;
        } finally {
            document.body.removeChild(textarea);
        }
    }
}

$('#navigation_bar').on('click', '#btn_sign_in', (e) => {
    // Reset the appearance in case there was a previous login error
    $('#user_email').css({ 'color':'black', 'font-weight':'normal' });
    $('#user_pass').css({ 'color':'black', 'font-weight':'normal' });

    e.preventDefault();

    const formData = $("#login_form").serializeArray();

    $.ajax({
        url : './cgi/login.cgi',
        type: "POST",
        data : formData,
        dataType:"json"
    }).done((data) => {
        // login actions here
        //  0 - The user name wasn't found at all.
        if (data.session_id == 0) {
            $('#user_email').focus();
            $('#user_email').css({ 'color':'red', 'font-weight':'bold' });

        // -1 - User was found, but the password was incorrect
        } else if (data.session_id == -1) {
            $('#user_pass').focus();
            $('#user_pass').css({ 'color':'red', 'font-weight':'bold' });
            $('#user_pass').parent().addClass('has-error');

        //  ? - Any other value is the session ID
        } else {
            CURRENT_USER.email = $('#user_email').val();
            CURRENT_USER.user_name = data.name;
            CURRENT_USER.session_id = data.session_id;
            CURRENT_USER.profile = data['gear_default_domain'];
            CURRENT_USER.id = data.user_id;

            CURRENT_USER.is_admin = data['is_admin'] == 1;

            $('span.user_logged_in').text(CURRENT_USER.user_name);
            Cookies.set('gear_session_id', CURRENT_USER.session_id, { expires: 7 });
            session_id = Cookies.get('gear_session_id');
            Cookies.set('gear_default_domain', CURRENT_USER.profile);

            // do we process the current page or reload?
            if (document.URL.includes("dataset_explorer.html")) {
                location.reload();
            } else {
                handle_login_ui_updates();
            }
        }
    }).fail((jqXHR, textStatus, errorThrown) => {
        alert(`Failure!  Status: (${textStatus}) Error: (${errorThrown})`);
    });
});

$('#navigation_bar').on('keyup', '#user_pass', function(e){
    //reset password input when user types
    $('#user_pass').parent().removeClass('has-error');
    $('#user_pass').css({ 'color':'black', 'font-weight':'normal' });
    $('.popover').popover('hide');
});

$('#navigation_bar').on('click', '#btn_sign_out', function(e){
    $('#user_email').val('');
    $('#user_pass').val('');
    $('#search_gene_symbol').focus();

    Cookies.remove('gear_session_id');
    // these should be cleared upon page refresh but adding for sanity's sake
    session_id = undefined;
    CURRENT_USER = undefined;
    e.preventDefault();

    // prevents popover from popping up just before changing url. bootstrap bug?
    $('#search_gene_symbol').popover('dispose');
    // just redirect to home page
    window.location.replace('./index.html');
});

function handle_login_ui_updates() {
    if (CURRENT_USER.session_id == null) {
        // these are the pages which require a login
        if (document.URL.includes("upload_dataset.html") ||
            document.URL.includes("analyze_dataset.html") ||
            document.URL.includes("projection.html") ||
            document.URL.includes("user_profile.html") ||
            document.URL.includes("upload_epigenetic_data.html") ||
            document.URL.includes("dataset_curator.html") ||
            document.URL.includes("gene_cart_manager.html") ||
            document.URL.includes("multigene_curator.html") ||
            document.URL.includes("epiviz_panel_designer.html")) {
            $('div#login_warning').show();
            $('div#login_checking').hide();
            $('div#main_content').hide();
        }

        if (document.URL.includes("manual.html")) {
            $('a#user_guide_link').parent().addClass('active');

        } else if (document.URL.includes("contact.html")) {
            $('a#comment_link').parent().addClass('active');

        } else if (document.URL.includes("projection.html")) {
            populate_dataset_selection();
        } else if (document.URL.includes("dataset_explorer.html")) {
            // these are defined in dataset_explorer.js
            load_preliminary_data();
            $('div#login_checking').hide();
            $("#controls_profile_c").remove();
        } else if (document.URL.includes("gene_cart_manager.html")) {
            $('div#login_checking').hide();
        }

        $('#login_controls').show();
        //TODO: This is ugly, but hide() doesn't work here.
        $('#loggedin_controls').attr("style", "display: none !important");
        return;
    }

    $('#loggedin_controls').show();
    //TODO: This is ugly, but hide() doesn't work here.
    $('#login_controls').attr("style", "display: none !important");

    if (document.URL.includes("compare_datasets.html")) {
        $("#create_gene_cart").prop("disabled", false);
        $("#create_gene_cart").attr("title", "Create a gene cart with these genes");

    } else if (document.URL.includes("contact.html")) {
        $('a#comment_link').parent().addClass('active');

    } else if (document.URL.includes("create_account.html")) {
        $('#account_creation_form_c').hide();
        $('#account_already_created_c').show();

    } else if (document.URL.includes("dataset_explorer.html")) {
        // these are defined in dataset_explorer.js
        load_preliminary_data();
        $("#controls_profile_nouser_c").remove();
        $("#your_dataset_filter").show();

    } else if (document.URL.includes("gene_cart_manager.html")) {
        // these are defined in dataset_explorer.js
        load_preliminary_data();

    } else if (document.URL.includes("manual.html")) {
        $('a#user_guide_link').parent().addClass('active');

    }

    if (document.URL.includes("upload_dataset.html") ||
        document.URL.includes("dataset_explorer.html") ||
        document.URL.includes("gene_cart_manager.html") ||
        document.URL.includes("analyze_dataset.html") ||
        document.URL.includes("projection.html") ||
        document.URL.includes("user_profile.html") ||
        document.URL.includes("upload_epigenetic_data.html") ||
        document.URL.includes("dataset_curator.html") ||
        document.URL.includes("multigene_curator.html") ||
        document.URL.includes("epiviz_panel_designer.html")) {
        $('div#login_warning').hide();
        $('div#login_checking').hide();
        $('div#main_content').show();
        $('input#session_id').val(CURRENT_USER.session_id);
    }
}

$(document).on('click', 'button#submit_forgot_pass_email', function(e) {
    e.preventDefault();
    $('p#forgot_pass_instruct, p#forgot_pass_warning, input#forgot_pass_email, button#submit_forgot_pass_email').hide();
    $('#submit_wait_c').show();

    const email = $('input#forgot_pass_email').val();

    //send a trimmed current url. Helps to be more dynamic for easier testing
    const current_url = window.location.href;
    const current_page = current_url.lastIndexOf("/");
    const destination_page = `${current_url.substring(0, current_page)}/index.html`;
    $.ajax({
        url: './cgi/send_email.cgi',
        type: 'POST',
        data: { 'email': email, 'scope': 'forgot_password', 'destination_page': destination_page },
        dataType: 'json',
        success(data) {
            $('#submit_wait_c').hide();
            if (data.success == 1) {
                $('#forgot_pass_c').hide();
                $('#forgot_pass_requested_c').show();
                return;
            }
            $('p#forgot_pass_instruct').hide();
            $('input#forgot_pass_email').val('');
            $('p#forgot_pass_warning, input#forgot_pass_email, button#submit_forgot_pass_email').show();
        },
        error(jqXHR, textStatus, errorThrown) {
            display_error_bar(`${jqXHR.status} ${errorThrown.name}`, 'Failure to submit forgotten password form');
        }
    });
});

//Prevent clicking on forgot pass button from going somewhere
$('#navigation_bar').on('click', '#forgot_password_link', function(e) {
    e.preventDefault();
    load_forgot_password();
});

function load_forgot_password() {
    $('#forgot_password_link[data-toggle="popover"]').popover({
        animation: true,
        trigger: 'click',
        title: "Forgot your password? <button type='button' id='dismiss_fp' class='close' style='margin:-2px;'>&times;</button>",
        html: true,
        placement: 'auto'
    }).on('show.bs.popover', function(e) {
        // Helped to insert the share_url more reliably: http://stackoverflow.com/a/25885326/2900840
        const el = $(e.target);

        // Trick bootstrap. Remove the title. bootstrap looks to next option for a title - that's above.
        $('#forgot_password_link').attr('data-original-title', "Forgot your password? <button type='button' id='dismiss_fp' class='close' style='margin:-2px;'>&times;</button>");

        // Help to close other popovers: http://stackoverflow.com/a/34320956/2900840
        el.data("bs.popover")._activeTrigger.click = false;

        const html = "<div id='forgot_pass_c' class='text-center'>" +
                "<p id='forgot_pass_instruct'>Enter your email and we will send you a link to reset your password.</p>" +
                "<p id='forgot_pass_warning' class='text-warning shown_later'>We could not find that email. Please check what you entered and try again.</p>" +
                "<input type='text' id='forgot_pass_email' class='form-control mb-1' value=''>" +
                "<button id='submit_forgot_pass_email' class='btn btn-success' data-dismiss='popover' value='submit'>Submit</button>" +
                "<div id='submit_wait_c' style='display:none;text-align:center;'>" +
                "<img style='height:60px;' src='img/loading_search.gif' alt='submitting'>" +
                "<h2>Please wait</h2>" +
                "<p>We're looking up your email address.</p>" +
                "</div>" +
                "</div>" +
                "<div id='forgot_pass_requested_c' class='text-center shown_later'>" +
                "<h1><span style='cursor:none;' class='glyphicon glyphicon-send'></span></h1>" +
                "<h2>Request sent.</h2>" +
                "<p>Please check your email.</p>" +
                "</div>";

        //insert content into share popover
        el.attr('data-content', html);
    });
}

// Close the Forgot Password popover
$(document).on('click', 'button#dismiss_fp', function() {
    $('.popover').popover('hide');

    // Add back the tooltip
    $('#forgot_password_link').attr('data-original-title', 'Reset your password');
});

$(document).on('click', 'button.click-once', function(e) {
    /*
       Generic listener to make sure buttons get disabled when clicked.  It is up
       to the secondary event to re-enable it.
    */
    $(this).attr("disabled", true);
});

// If a video link is clicked, open the Youtube video in a new window
$(document).on("click", ".js-video-link", (e) => {

    if (document.URL.includes("dataset_explorer.html")) {
        window.open("https://youtu.be/jZVCF2Yqm4M");
    }

    if (document.URL.includes("compare_datasets.html")) {
        window.open("https://youtu.be/msmppWq6XrQ");
    }

    if ((document.URL.includes("index.html")) || (['/', 'index.html'].includes(location.pathname))) {
        window.open("https://youtu.be/sr_kvm7W4OE");
    }

    return false;
});

// If user is on one of these pages, show the video walkthrough link.  Else hide it (since we do not have a video for it).
function show_video_link() {
    if ((document.URL.includes("dataset_explorer.html"))
        || (document.URL.includes("compare_datasets.html"))
        || (document.URL.includes("index.html"))
        || (['/', 'index.html'].includes(location.pathname))) {
            return true;
        }

    return false;
}

// If user guide link is clicked, navigate to manual page and potentially to the anchor
$(document).on("click", ".js-user-guide-link", (e) => {
    let anchor;
    let manual_url = "/manual.html"

    if (document.URL.includes("upload_dataset.html")) {anchor="uploading"}
    if (document.URL.includes("analyze_dataset.html")) {anchor="workbench"}
    if (document.URL.includes("upload_epigenetic_data.html")) {anchor="epiviz"}
    if (document.URL.includes("dataset_curator.html")) {anchor="curation"}
    if (document.URL.includes("dataset_explorer.html")) {anchor="profiles"}
    //if (document.URL.includes("gene_cart_manager.html")) {}
    if (document.URL.includes("compare_datasets.html")) {anchor="compare"}
    if (document.URL.includes("multigene_curator.html")) {anchor="mg-curation"}
    if (document.URL.includes("epiviz_panel_designer.html")) {anchor="epiviz"}

    if (anchor) {
        manual_url += `?doc=${anchor}`
    }

    window.open(manual_url);
    return false;
});

// This function takes the contents of an HTML table and formats it as a tab
//  delimited string with rows and then offers it as a file download, with
//  excel file extension.
// Expects the table to have a proper thead with tr/th elements, and tbody with tr/td
function download_table_as_excel(table_id, filename) {
    table_str = '';

    $('#' + table_id + ' thead tr th').each(function() {
        table_str += $(this).text() + "\t";
    });
    table_str = table_str.trim() + "\n";

    $('#' + table_id + ' tbody tr').each(function() {
        var rows = $(this).find('td');

        rows.each(function() {
            table_str += $(this).text() + "\t";
        });

        table_str = table_str.trim() + "\n";
    });

    var element = document.createElement('a');
    element.setAttribute('href', 'data:text/plain;charset=utf-8,' + encodeURIComponent(table_str));
    element.setAttribute('download', filename);
    element.style.display = 'none';
    document.body.appendChild(element);
    element.click();
    document.body.removeChild(element);
}

// Generates an RFC4122 version 4 compliant UUID
function uuid() {
  return 'xxxxxxxx-xxxx-xxxx-yxxx-xxxxxxxxxxxx'.replace(/[xy]/g, function(c) {
    const r = Math.random() * 16 | 0, v = c == 'x' ? r : (r & 0x3 | 0x8);
    return v.toString(16);
  });
}


// From: http://stackoverflow.com/a/21903119/1368079
// Use example:
//   if URL is: http://dummy.com/?technology=jquery&blog=jquerybyexample
//   then:      var tech = getUrlParameter('technology');
//              var blog = getUrlParameter('blog');
const getUrlParameter = function getUrlParameter(sParam) {
    const sPageURL = decodeURIComponent(window.location.search.substring(1));
    const sURLVariables = sPageURL.split('&');
    let sParameterName;
    let i;

    for (i = 0; i < sURLVariables.length; i++) {
        sParameterName = sURLVariables[i].split('=');

        if (sParameterName[0] === sParam) {
            return sParameterName[1] === undefined ? true : sParameterName[1];
        }
    }
};

// error should be html message for user. Example: error = '<p>You cannot do that.</p>'
function display_error_bar(msg, detail) {
    $('.alert-container').html('<div class="alert alert-danger alert-dismissible" role="alert">' +
      '<button type="button" class="close close-alert" data-dismiss="alert" aria-label="Close"><span aria-hidden="true">&times;</span></button>' +
      '<p class="alert-message">' +
      '<strong>Fail. </strong> Sorry, something went wrong: ' + detail + '  Please contact us with this message for help.' +
      '</p>' +
      '<p style="text-align: right;">(<em>Error: ' + msg + '</em>)</p>' +
      '</div>').show();
}

/* https://naveensnayak.com/2013/06/26/dynamically-loading-css-and-js-files-using-jquery/ */
loadCSS = (href) => {
    const cssLink = $(`<link rel='stylesheet' type='text/css' href='${href}'>`);
    $("head").append(cssLink);
};

function image_loaded(img) {
    if (!img.complete) {
        return false;
    }

    if (typeof img.naturalWidth != "undefined" && img.naturalWidth == 0) {
        return false;
    }

    return true;
}
