var search_results = [];
var layouts = []; //populated on load by loads_layouts();
var gene_carts = [];

// keyed first on dataset #, then body site
// like site_nodes['dataset1'][5] = array of path nodes
var site_nodes = {};

// key dataset_id, value = Snap paths
var svgs = {};

var SCROLLBAR_DRAWN = false;
var GO_TERM_SCROLLBAR_DRAWN = false;
var AT_FIRST_MATCH_RECORD = false;
var AT_LAST_MATCH_RECORD  = false;
var PREVIOUS_SELECTED_RECORD_NUM = null;
var SCORING_METHOD = 'gene';
var SELECTED_GENE = null;

var COLORS = {};
var SELECTED_COLOR_SET = 'linear_purple';
COLORS['netherlands'] = ['rgb(0,0,255)','rgb(2,2,255)','rgb(4,4,255)','rgb(6,6,255)','rgb(8,8,255)','rgb(10,10,255)','rgb(12,12,255)','rgb(14,14,255)','rgb(16,16,255)','rgb(18,18,255)','rgb(20,20,255)','rgb(22,22,255)','rgb(24,24,255)','rgb(26,26,255)','rgb(28,28,255)','rgb(30,30,255)','rgb(32,32,255)','rgb(34,34,255)','rgb(36,36,255)','rgb(38,38,255)','rgb(40,40,255)','rgb(42,42,255)','rgb(44,44,255)','rgb(46,46,255)','rgb(48,48,255)','rgb(50,50,255)','rgb(52,52,255)','rgb(54,54,255)','rgb(56,56,255)','rgb(58,58,255)','rgb(60,60,255)','rgb(62,62,255)','rgb(64,64,255)','rgb(66,66,255)','rgb(68,68,255)','rgb(70,70,255)','rgb(72,72,255)','rgb(74,74,255)','rgb(76,76,255)','rgb(78,78,255)','rgb(80,80,255)','rgb(82,82,255)','rgb(84,84,255)','rgb(86,86,255)','rgb(88,88,255)','rgb(90,90,255)','rgb(92,92,255)','rgb(94,94,255)','rgb(96,96,255)','rgb(98,98,255)','rgb(100,100,255)','rgb(102,102,255)','rgb(104,104,255)','rgb(106,106,255)','rgb(108,108,255)','rgb(110,110,255)','rgb(112,112,255)','rgb(114,114,255)','rgb(116,116,255)','rgb(118,118,255)','rgb(120,120,255)','rgb(122,122,255)','rgb(124,124,255)','rgb(126,126,255)','rgb(128,128,255)','rgb(130,130,255)','rgb(132,132,255)','rgb(134,134,255)','rgb(136,136,255)','rgb(138,138,255)','rgb(140,140,255)','rgb(142,142,255)','rgb(144,144,255)','rgb(146,146,255)','rgb(148,148,255)','rgb(150,150,255)','rgb(152,152,255)','rgb(154,154,255)','rgb(156,156,255)','rgb(158,158,255)','rgb(160,160,255)','rgb(162,162,255)','rgb(164,164,255)','rgb(166,166,255)','rgb(168,168,255)','rgb(170,170,255)','rgb(172,172,255)','rgb(174,174,255)','rgb(176,176,255)','rgb(178,178,255)','rgb(180,180,255)','rgb(182,182,255)','rgb(184,184,255)','rgb(186,186,255)','rgb(188,188,255)','rgb(190,190,255)','rgb(192,192,255)','rgb(194,194,255)','rgb(196,196,255)','rgb(198,198,255)','rgb(200,200,255)','rgb(202,202,255)','rgb(204,204,255)','rgb(206,206,255)','rgb(208,208,255)','rgb(210,210,255)','rgb(212,212,255)','rgb(214,214,255)','rgb(216,216,255)','rgb(218,218,255)','rgb(220,220,255)','rgb(222,222,255)','rgb(224,224,255)','rgb(226,226,255)','rgb(228,228,255)','rgb(230,230,255)','rgb(232,232,255)','rgb(234,234,255)','rgb(236,236,255)','rgb(238,238,255)','rgb(240,240,255)','rgb(242,242,255)','rgb(244,244,255)','rgb(246,246,255)','rgb(248,248,255)','rgb(250,250,255)','rgb(252,252,255)','rgb(255,253,253)','rgb(255,251,251)','rgb(255,249,249)','rgb(255,247,247)','rgb(255,245,245)','rgb(255,243,243)','rgb(255,241,241)','rgb(255,239,239)','rgb(255,237,237)','rgb(255,235,235)','rgb(255,233,233)','rgb(255,231,231)','rgb(255,229,229)','rgb(255,227,227)','rgb(255,225,225)','rgb(255,223,223)','rgb(255,221,221)','rgb(255,219,219)','rgb(255,217,217)','rgb(255,215,215)','rgb(255,213,213)','rgb(255,211,211)','rgb(255,209,209)','rgb(255,207,207)','rgb(255,205,205)','rgb(255,203,203)','rgb(255,201,201)','rgb(255,199,199)','rgb(255,197,197)','rgb(255,195,195)','rgb(255,193,193)','rgb(255,191,191)','rgb(255,189,189)','rgb(255,187,187)','rgb(255,185,185)','rgb(255,183,183)','rgb(255,181,181)','rgb(255,179,179)','rgb(255,177,177)','rgb(255,175,175)','rgb(255,173,173)','rgb(255,171,171)','rgb(255,169,169)','rgb(255,167,167)','rgb(255,165,165)','rgb(255,163,163)','rgb(255,161,161)','rgb(255,159,159)','rgb(255,157,157)','rgb(255,155,155)','rgb(255,153,153)','rgb(255,151,151)','rgb(255,149,149)','rgb(255,147,147)','rgb(255,145,145)','rgb(255,143,143)','rgb(255,141,141)','rgb(255,139,139)','rgb(255,137,137)','rgb(255,135,135)','rgb(255,133,133)','rgb(255,131,131)','rgb(255,129,129)','rgb(255,127,127)','rgb(255,125,125)','rgb(255,123,123)','rgb(255,121,121)','rgb(255,119,119)','rgb(255,117,117)','rgb(255,115,115)','rgb(255,113,113)','rgb(255,111,111)','rgb(255,109,109)','rgb(255,107,107)','rgb(255,105,105)','rgb(255,103,103)','rgb(255,101,101)','rgb(255,99,99)','rgb(255,97,97)','rgb(255,95,95)','rgb(255,93,93)','rgb(255,91,91)','rgb(255,89,89)','rgb(255,87,87)','rgb(255,85,85)','rgb(255,83,83)','rgb(255,81,81)','rgb(255,79,79)','rgb(255,77,77)','rgb(255,75,75)','rgb(255,73,73)','rgb(255,71,71)','rgb(255,69,69)','rgb(255,67,67)','rgb(255,65,65)','rgb(255,63,63)','rgb(255,61,61)','rgb(255,59,59)','rgb(255,57,57)','rgb(255,55,55)','rgb(255,53,53)','rgb(255,51,51)','rgb(255,49,49)','rgb(255,47,47)','rgb(255,45,45)','rgb(255,43,43)','rgb(255,41,41)','rgb(255,39,39)','rgb(255,37,37)','rgb(255,35,35)','rgb(255,33,33)','rgb(255,31,31)','rgb(255,29,29)','rgb(255,27,27)','rgb(255,25,25)','rgb(255,23,23)','rgb(255,21,21)','rgb(255,19,19)','rgb(255,17,17)','rgb(255,15,15)','rgb(255,13,13)','rgb(255,11,11)','rgb(255,9,9)','rgb(255,7,7)','rgb(255,5,5)','rgb(255,3,3)','rgb(255,1,1)','rgb(255,0,0)'];
COLORS['linear_purple'] = ['#f9f9f9','#f8f8f8','#f7f7f6','#f6f6f5','#f5f5f4','#f3f5f2','#f1f4f2','#eff3f2','#eef2f2','#eceff2','#ebebf1','#ebe9f1','#ede7f0','#f0e6ef','#efe4ea','#efe2e4','#efe2e4','#efe2e4','#eee1e3','#eedfe1','#eddee0','#ecdcdf','#ecdbde','#ebdadc','#ead8db','#ead7da','#e9d6d9','#e8d4d7','#e8d3d6','#e7d1d5','#e6d0d4','#e6ced2','#e5cdd1','#e5cdd1','#e5cdd1','#e4ccd0','#e4cacf','#e3c9ce','#e2c7cd','#e2c6cc','#e1c5cb','#e0c3ca','#e0c2c9','#dfc0c8','#debfc7','#debec6','#ddbcc5','#dcbbc4','#dcb9c3','#dbb8c2','#dbb8c2','#dbb8c2','#dab7c1','#dab5c0','#d9b4bf','#d8b3be','#d8b1bd','#d7b0bd','#d6afbc','#d6adbb','#d5acba','#d4abb9','#d4a9b8','#d3a8b7','#d2a7b7','#d2a5b6','#d1a4b5','#d1a4b5','#d1a4b5','#d0a3b4','#d0a1b4','#cfa0b3','#ce9fb2','#ce9eb1','#cd9cb1','#cc9bb0','#cc9aaf','#cb99af','#ca97ae','#ca96ad','#c995ad','#c894ac','#c892ac','#c791ab','#c791ab','#c791ab','#c690aa','#c68faa','#c58ea9','#c48ca8','#c48ba8','#c38aa7','#c289a7','#c288a6','#c187a5','#c086a5','#c085a4','#bf83a4','#be82a3','#be81a3','#bd80a2','#bd80a2','#bd80a2','#bc7fa2','#bc7ea1','#bb7da1','#ba7ba0','#ba7aa0','#b9799f','#b8789f','#b8779f','#b7769e','#b6759e','#b6749d','#b5729d','#b4719d','#b4709c','#b36f9c','#b36f9c','#b36f9c','#b26e9c','#b26d9b','#b16c9b','#b06b9b','#af6a9a','#af699a','#ae689a','#ad6799','#ac6699','#ac6599','#ab6498','#aa6398','#a96298','#a96197','#a86097','#a86097','#a86097','#a85f97','#a75d97','#a75c97','#a75b97','#a65997','#a65896','#a55896','#a45796','#a35695','#a25595','#a15494','#a15394','#a05394','#9f5293','#9e5193','#9e5193','#9e5193','#9d5093','#9d4f93','#9c4e92','#9b4d92','#9b4d92','#9a4c92','#994b92','#994a91','#984991','#974891','#974791','#964791','#954690','#954590','#944490','#944490','#944490','#934390','#934290','#92428f','#91418f','#91408f','#903f8f','#8f3e8f','#8f3d8e','#8e3d8e','#8d3c8d','#8c3b8d','#8b3a8c','#893a8b','#88398b','#87388a','#87388a','#87388a','#863789','#853689','#843688','#823587','#813487','#803386','#7f3385','#7e3285','#7c3184','#7b3183','#7a3083','#792f82','#772e81','#762e81','#752d80','#752d80','#752d80','#742c7f','#732c7f','#712b7e','#702a7d','#6f2a7d','#6e297c','#6d287b','#6b287b','#6a277a','#692679','#682679','#672578','#652477','#642477','#632376','#632376','#632376','#622275','#612275','#5f2174','#5e2073','#5d2073','#5c1f72','#5b1f71','#591e71','#581d70','#571d6f','#561c6f','#551c6e','#531b6d','#521b6d','#511a6c','#511a6c','#511a6c','#50196b','#4f196b','#4e196a','#4c1869','#4b1869','#4a1768','#491767','#481667','#471666','#461565','#441565','#431464','#421463','#411363','#401362','#401362'];
COLORS['custom_relative'] = [];
COLORS['custom_absolute'] = [];
var COLOR_MODE = 'default';

var share_id = null; //from permalink
var permalinked_dataset_id = null; //holds dataset_id obtained from load_dataset_frames()

var annotation_panel = new FunctionalAnnotationPanel();
var dataset_collection_panel = new DatasetCollectionPanel();

var search_result_postselection_functions = [];

window.onload=function() {
    // check if the user is already logged in
    check_for_login();

    $('#highlighted_dataset').load('./include/by_domain/' + SITE_PREFS['domain_label'] + '/index_highlighted_dataset.html');

    // Was a permalink found?
    var share_id = getUrlParameter('share_id');
    if (share_id) {
        //hide site_into and display the permalink message
        $('#intro_content').hide();
        $('#viewport_intro').children().hide();
        $('#searching_indicator_c').hide();

        $('#leftbar_main').show();
        $('#permalink_intro_c').show();
        // validate the share_id. runs load_dataset_frames() on success
        validate_permalink(share_id, 'permalink');
    } else {
        get_index_info();
        load_layouts();
        load_gene_carts();
    }

    // Was help_id found?
    var help_id = getUrlParameter('help_id');
    if (help_id) {
        validate_help_id(help_id);
    }

    var permalinked_gene_symbol = getUrlParameter('gene_symbol');
    var permalinked_gsem = getUrlParameter('gene_symbol_exact_match');
    if (permalinked_gene_symbol) {
        $("#search_gene_symbol_intro").val(permalinked_gene_symbol);

        if (permalinked_gsem) {
            set_exact_match('on');
        }

        sleep(1000).then(() => {
            $('#intro_search_icon').trigger('click');
            // clear any open tooltips
            $('[data-toggle="tooltip"], .tooltip').tooltip("hide");
        })
    }

    // The search button starts out disabled, make sure it gets re-enabled.
    $("button#submit_search").prop( "disabled", false );

    $('#exact_match_icon').click(function() {
        // handle if it was already in the on position
        if ( $('#exact_match_input').prop("checked") ) {
            set_exact_match('off');
        } else {
            // else it was off, so turn it on
            set_exact_match('on');
        }
    });

    $('#intro_search_form').on('submit', function(e) {
        // TODO: It makes sense to remove/destroy those elements we aren't showing after a search
        e.preventDefault();
        $("#search_gene_symbol").val( $("#search_gene_symbol_intro").val() );
        $('#intro_content').hide();

        $("#leftbar_main").show();
        $("#viewport_main").show();

        // Set to 1 to add anchor to page, treating as a new page.
        // Useful for retrieving original home page layout if user hits 'Back' button.
        // Source - https://www.webdesignerdepot.com/2013/03/how-to-manage-the-back-button-with-javascript/
        location.hash = 1;

        // fire the true search button
        $("#submit_search").trigger( "click" );
    });

    $('#intro_search_icon').click(function() {
        $('#intro_search_form').submit();
    });

    $('#launcher_comparison_tool').click(function() {
        window.location.replace('./compare_datasets.html');
    });

    $('#launcher_dataset_manager').click(function() {
        window.location.replace('./dataset_manager.html');
    });

    $('#launcher_manual').click(function() {
        window.location.replace('./manual.html');
    });

    $('#launcher_uploader').click(function() {
        window.location.replace('./upload_dataset.html');
    });

    $('#launcher_workbench').click(function() {
        window.location.replace('./analyze_dataset.html');
    });

    // add post-page load listeners
    $( "#dataset_zoomed_zoom_out_control" ).click(function() {
        zoom_out_dataset();
    });

    $(document).on('click', '.domain_choice_c', function() {
        dataset_collection_panel.set_layout($(this).data('profile-id'), $(this).data('profile-label'), true);
    });

    $( document ).on("click", ".scope_choice", function() {
        SCORING_METHOD = $(this).data('choice');
        if (SELECTED_GENE !== null) {
            select_search_result($(SELECTED_GENE));
        }
    });

    // track the mouse movement so we can display scoring tooltips
    $( document ).on( "mousemove", function( event ) {
        // Positioning for dataset_grid tips
        // Why is this pixel adjustment necessary?
        xpos = event.pageX - 240;
        ypos = event.pageY - 130 - 30;
        $("#tip").css("left", xpos + "px" );
        $("#tip").css("top" , ypos + "px" );
    });
};

// Load original home page if back button was clicked from a gene search results page.
window.onhashchange = function() {
    // On original index page
    if (location.hash.length === 0) {
        SELECTED_GENE = null;
        $('#intro_content').show();
        $('#intro_content_c').show();

        $("#leftbar_main").hide();
        $("#viewport_main").hide();
    }
    // Navigating to gene search result version of index page makes url look like 'localhost/#'
}

function get_index_info() {
    $.ajax({
        url: './cgi/get_index_info.cgi',
        type: 'GET',
        dataType: 'json',
        success: function(data, textStatus, jqXHR) {
            /* Not currentlyd oing this
            var tmpl = $.templates("#tmpl_newest_datasets");
            var thtml = tmpl.render(data['newest_public_datasets']['datasets']);
            $("#new_datasets").html(thtml);
            */

            $('#stats_dataset_count').text(data['dataset_count'])
            $('#stats_user_count').text(data['user_count'])
        },
        error: function(jqXHR, textStatus, errorThrown) {
            display_error_bar(jqXHR.status + ' ' + errorThrown.name, 'Error getting index info.');
        }
    });
};

//Check help_id is valid. For Forgotten Password
function validate_help_id(help_id) {
    $.ajax({
        url: './cgi/validate_help_id.cgi',
        type: 'POST',
        data: {'help_id': help_id},
        dataType: 'json',
        success: function(data, textStatus, jqXHR) {
            if ( data['success'] == 1 ) {
                // Add help_id to form
                $('#user_help_id').val(help_id);

                // Greet user by name (a subtle confirmation they know it's their account)
                if (data['user_name'].length > 0) {
                    var user_first_name = ' ' + data['user_name'].split(' ')[0] + '!';
                    $('#forgot_password_user_name').text(user_first_name);
                }

                $('#forgot_password_modal').modal('show');
            } else {
                // Invalid help_id, display invalid message
                $("#valid_forgot_pass_modal_body_c").hide();
                $("#save_user_new_pass").hide();
                $("#invalid_forgot_pass_modal_body_c").show();

                $('#forgot_password_modal').modal('show');
            }
        },
        error: function(jqXHR, textStatus, errorThrown) {
            display_error_bar(jqXHR.status + ' ' + errorThrown.name, 'Error validating help ID');
        }
    });
};

// Disable 2nd password input until 1st is populated
$(document).on('keydown', 'input#user_new_pass_1', function(){
    if ( $(this).val().length > 1 ) {
        $("input#user_new_pass_2").prop("disabled", false);
    } else {
        $("input#user_new_pass_2").prop("disabled", true);
    }
});

// Disable Save password button until 1st and 2nd inputs match
$(document).on('keyup', 'input#user_new_pass_2', function(){
    var pass_1 = $('input#user_new_pass_1').val();

    if ( $('input#user_new_pass_2').val() == pass_1 ) {
        $('button#save_user_new_pass').prop('disabled', false);
    } else {
        $('button#save_user_new_pass').prop('disabled', true);
    }
});

// Submit new password
$(document).on('click', 'button#save_user_new_pass', function(){
    // Hide password form and show waiting
    $('#valid_forgot_pass_modal_body_c').hide();
    $('#forgot_pass_modal_footer').hide();
    $('#saving_forgot_pass_modal_body_c').show();

    var help_id = $('input#user_help_id').val();
    var new_password = $('input#user_new_pass_2').val();
    $.ajax({
        url: './cgi/save_user_account_changes.cgi',
        type: 'POST',
        data: { 'help_id': help_id, 'new_password': new_password, 'scope': 'password'},
        dataType: 'json',
        success: function(data, textStatus, jqXHR) {
            if (data['success'] == 1) {
                // Hide waiting and show success
                $('#saving_forgot_pass_modal_body_c').hide();
                $('#success_forgot_pass_modal_body_c').show();

                // Redirect to home page
                setInterval(function(){
                    window.location.replace('./index.html');
                }, 2000);
            } else {
                $('.alert-container').html('<div class="alert alert-danger alert-dismissible" role="alert">' +
                    '<button type="button" class="close close-alert" data-dismiss="alert" aria-label="Close"><span aria-hidden="true">&times;</span></button>' +
                    '<p class="alert-message"><strong>Oops! </strong> ' + data["error"] + '</p></div>').show();
            }
        },
        error: function (jqXHR, textStatus, errorThrown) {
            display_error_bar(jqXHR.status + ' ' + errorThrown.name, 'Error saving new user password');
        }
    });//end ajax
});

function validate_permalink(share_id, scope) {
    $.ajax({
        url : './cgi/validate_share_id.cgi',
        type: "POST",
        data : { 'share_id': share_id, 'scope': scope },
        dataType:"json",
        success: function(data, textStatus, jqXHR) {
            if ( data['success'] == 1 ) {
                // query the db and load the images, including permalink dataset
                dataset_collection_panel.load_frames({ share_id });

            } else {
                // query the db and load the images
                dataset_collection_panel.load_frames();
                $('.alert-container').html('<div class="alert alert-danger alert-dismissible" role="alert">' +
                    '<button type="button" class="close close-alert" data-dismiss="alert" aria-label="Close"><span aria-hidden="true">&times;</span></button>' +
                    '<p class="alert-message"><strong>Oops! </strong> ' + data["error"] + '</p></div>').show();
            }
        },
        error: function (jqXHR, textStatus, errorThrown) {
            display_error_bar(jqXHR.status + ' ' + errorThrown.name, 'Error validating share ID');
        }
    });
}

function load_layouts() {
    var d = new $.Deferred();
    var session_id = Cookies.get('gear_session_id');
    var layout_share_id = getUrlParameter('layout_id');

    //organize user and domain profiles in x-editable format
    $.ajax({
        url: './cgi/get_user_layouts.cgi',
        type: 'post',
        data: { 'session_id': session_id, 'layout_share_id': layout_share_id },
        dataType: 'json',
        success: function(data, textStatus, jqXHR) {
            /*
              Priority of displayed profile:
                0.  Passed layout ID via layout_id URL parameter
                1.  Cookie value
                2.  User's DB-saved value (when they go to a machine, and there's no cookie)
                3.  Admin's active domain
             */

            var domain_profiles = [];
            var user_profiles = [];

            var active_layout_id = null;
            var active_layout_label = null;

            // Pass through once to sort domains from user profiles AND see if it matches a shared layout
            $.each(data['layouts'], function(i, item){
                if ( item['is_domain'] == 1 ) {
                    domain_profiles.push({value: item['id'], text: item['label'] });
                } else {
                    user_profiles.push({value: item['id'], text: item['label'] });
                }

                if (item['share_id'] == layout_share_id) {
                    active_layout_id = item['id'];
                    active_layout_label = item['label'];
                }
            });

            // pass through again and look for one set by a cookie
            if (active_layout_id == null) {
                $.each(data['layouts'], function(i, item) {
                    if (item['label'] == CURRENT_USER.profile) {
                        active_layout_id = item['id'];
                        active_layout_label = item['label'];
                        return false;
                    }
                });
            }

            // pass through again and look for one set as current by the user
            if (active_layout_id == null) {
                $.each(data['layouts'], function(i, item) {
                    if ( item['is_domain'] == 0 && item['is_current'] == 1 ) {
                        active_layout_id = item['id'];
                        active_layout_label = item['label'];
                        return false;
                    }
                });
            }

            // pass through again if no active layout was found for user and choose the admin's
            if (active_layout_id == null) {
                $.each(data['layouts'], function(i, item) {
                    if ( item['is_domain'] == 1 && item['is_current'] == 1 ) {
                        active_layout_id = item['id'];
                        active_layout_label = item['label'];
                        return false;
                    }
                });
            }

            var formattedData = [{text: 'Public profiles', children: domain_profiles}];

            if (user_profiles.length > 0) {
                formattedData.push({text: 'Your profiles', children: user_profiles});
            }

            var layout_items_tmpl = $.templates("#tmpl_intro_search_layouts");

            var domain_items_html = layout_items_tmpl.render(domain_profiles);
            $(domain_items_html).insertAfter('#intro_domains_header');

            var profile_items_html = layout_items_tmpl.render(user_profiles);
            $(profile_items_html).insertAfter('#intro_profiles_header');

            //Serves as source for #selected_profile editable
            layouts = formattedData;

            dataset_collection_panel.set_layout(active_layout_id, active_layout_label, true);

            d.resolve();
        },
        error: function (jqXHR, textStatus, errorThrown) {
            display_error_bar(jqXHR.status + ' ' + errorThrown.name, 'Error loading layouts.');
            d.fail();
        }
    });

    d.promise();
}

function load_gene_carts() {
  var d = new $.Deferred();
  var session_id = Cookies.get('gear_session_id');

  if (!session_id) {
      //User is not logged in. Hide gene carts container
      $("#selected_gene_cart_c").hide();
      d.resolve();
  } else {
      $("#selected_gene_cart_c").show(); //Show if hidden
      $.ajax({
        url: './cgi/get_user_gene_carts.cgi',
        type: 'post',
        data: { 'session_id': session_id },
        dataType: 'json',
        success: function(data, textStatus, jqXHR){ //source https://stackoverflow.com/a/20915207/2900840
            var user_gene_carts = [];
            var formattedData = [];

            if (data['gene_carts'].length > 0) {
                //User has some profiles
                $.each(data['gene_carts'], function(i, item){
                    user_gene_carts.push({value: item['id'], text: item['label'] });

                });

                formattedData = [{text: 'My gene carts', children: user_gene_carts }];
            } else {
                $("#selected_gene_cart_c").hide();
            }

            //Serves as source for #selected_gene_cart editable
            gene_carts = formattedData;
            d.resolve();
        },
        error: function (jqXHR, textStatus, errorThrown) {
            display_error_bar(jqXHR.status + ' ' + errorThrown.name);
            d.fail();
        }
      });
  }
  d.promise();
}

//Hide option menu when scope is changed.
$(document).on('click', '.scope_choice', function(){
    $('#toggle_options').popover('hide');
});

function populate_search_result_list(data) {
    // so we can display in sorted order.  javascript sucks like that.
    sorted_gene_syms = [];

    for (var key in data) {
        if (data.hasOwnProperty(key)) {
            sorted_gene_syms.push(key);
        }
    }

    sorted_gene_syms.sort();
    sorted_gene_syms_len = sorted_gene_syms.length

    var items = [];

    for (i = 0; i < sorted_gene_syms_len; i++) {
        gene_symbol = sorted_gene_syms[i];

        // Build search result html
        var gene_result_html = '<a class="list-group-item" data-gene_symbol="' + gene_symbol + '" href="#">' + gene_symbol;

        gene_result_html += '</a>';
        items.push(gene_result_html);
    }

    if (items.length == 0) {
        $('#search_results').text('No results found');
        $('#search_result_count').text(items.length);
    } else {
        $('#search_results').append( items.join('') );

        // the value here needs to match the max in gene_search.cgi
        if (items.length == 100) {
            $('#search_result_count').text('max:' + items.length);
        } else {
            $('#search_result_count').text(items.length);
        }
    }
}

var lastCall = 0;
function select_search_result(elm) {
    //TODO Prevents this function from being double-called by #gene_search_form.submit()
    var callTime = new Date().getTime();
    if (callTime - lastCall <= 500) {
      return false;
    }
    lastCall = callTime;

    SELECTED_GENE = $(elm);
    gene_sym = $(elm).data("gene_symbol");

    // remove coloring from other result links
    $('.list-group-item-active').removeClass('list-group-item-active');
    $(elm).addClass('list-group-item-active');

    var current_record_number = 0
    var match_record_number = null

    annotation_panel.annotation = search_results[gene_sym];
    annotation_panel.display_first_organism();

    // hide the intro, show the search result box
    if( $('#site_intro_c').is(':visible') ) {
        $('#site_intro_c').hide({easing: 'fade', duration: 400});
        $('#recent_updates_c').hide({easing: 'fade', duration: 400, complete: show_search_result_info_box});
    }

    dataset_collection_panel.update_by_search_result(search_results[gene_sym]);

    // call any plugin functions
    search_result_postselection_functions.forEach(function(f) {f()})
}

function isNumeric(n) {
  return !isNaN(parseFloat(n)) && isFinite(n);
}

function draw_graphs(expression_data, dataset_id) {
    //Gather expression data into lists of subplots
    var expressions_list = [];
    if (expression_data.length > 0) {
        var subplot_total = expression_data[0]['plots_total'];
        for (var i=0; i<subplot_total; i++){
            expressions_list[i] = [];
        }
        $.each(expression_data, function(i,item) {
            expressions_list[ item['plot_number']-1 ].push(item);
        });
    }

    var graph_type = null;
    // Here you might need to destroy any bar chart before redrawing
    // https://www.zingchart.com/docs/api/api-methods/#zingchart__exec__api__destroy
    var graph_container_target_id = null;

    //Get container for Bargraph
    if (dataset_id in bargraphs) {
        graph_container_target_id = "dataset_" + dataset_id + "_bargraph_standard";
        graph_type = 'bar';

        // Organize bargraph expression by graph position
        //    (older datasets lack 'graph_position', so should not be sorted)
        if (expressions_list.length > 0){
            $.each(expressions_list, function(i, list){
                if (list[0]['graph_position'] != null) {
                    list.sort(function(a, b) {
                        if (a.graph_position < b.graph_position) {return -1;}
                        if (a.graph_position > b.graph_position) {return 1;}
                        return 0;
                    });
                }
            });
        }
    }

    //Get container for Linegraph
    if (dataset_id in linegraphs) {
        graph_container_target_id = "dataset_" + dataset_id + "_linegraph_standard";
        graph_type = 'line';
    }

    //Destroy preexisting zingcharts
    if ( bargraphs[dataset_id] != null || linegraphs[dataset_id] != null ) {
        zingchart.exec(graph_container_target_id, 'destroy');
    }
    //get dataset's math_format
    var math_onload = $('#dataset_' + dataset_id + '_svg_c, #dataset_' + dataset_id + '_bargraph_standard').data('math-onload');

    var graphset_list = [];
    for (var plot=0, len=expressions_list.length; plot<len; plot++ ) {
        //Get subplot's expression values
        var expression = expressions_list[plot];

        var group_names = [];
        var groups = {};
        var max_series_size = 0;
        var series = [];
        var plot_title = null; //graphset var
        var plot_titles = [];
        var y_axis_label = null; // Will be the label for a single plot dataset (line graph only)
        var y_axis_labels = [];  // Will hold the labels for multi plot datasets
        for (var node=0, ex_len=expression.length; node<ex_len; node++) {

            // Bargraphs - Set group names and sec_labels
            if (graph_type == 'bar') {
                group_name = expression[node]['group_name'].replace(/_/g, ' ');
                sec_label = expression[node]['sec_label'];
            }

            // Linegraphs - Set group names. Sec_labels become Y-axis labels
            if (graph_type == 'line') {
                group_name = expression[node]['sec_label'].replace(/_/g, ' ');
                sec_label = null;
                if ($.inArray(expression[node]['group_name'], y_axis_labels)  == -1) {
                    y_axis_labels.push(expression[node]['group_name']);
                }
                y_axis_label = expression[node]['group_name'];
            }

            if (expression[node]['plot_title'] != null && expression[node]['plot_title'] != "_") {
                plot_title = expression[node]['plot_title'].replace(/_/g, ' ');
                plot_titles.push(plot_title);
            }

            if ($.inArray(group_name, group_names) == -1) {
                group_names.push(group_name);
                groups[group_name] = {values: [], errors: [], p_values: [], "data-xlabel": [], styles: []};
            }

            // graph is being redrawn (via math transformation)
            if ( expression[node]['new_val']) {
                expression_val = expression[node]['new_val'];
                std_dev_val = expression[node]['new_std_dev'];
            }
            // graph is being drawn (via selected gene)
            else {
                if ( math_onload == 'log2' ) {
                    expression_val = parseFloat(Math.log2(expression[node]['raw_val']).toFixed(4));
                    std_dev_val = null;
                    p_value = null;
                } else if ( math_onload == 'log10' ) {
                    expression_val = parseFloat(Math.log10(expression[node]['raw_val']).toFixed(4));
                    std_dev_val = null;
                    p_value = null;
                } else {
                    expression_val = expression[node]['raw_val'];
                    std_dev_val = expression[node]['std_dev'];
                    p_value = expression[node]['p_value'];
                }
            }

            groups[group_name]['values'].push(expression_val);
            groups[group_name]['errors'].push([std_dev_val]);
            groups[group_name]['data-xlabel'].push(sec_label);
            groups[group_name]['p_values'].push(p_value);

            if (groups[group_name]['values'].length > max_series_size) {
                max_series_size = groups[group_name]['values'].length;
                series.push(
                    { values: [], errors: [], p_values: [], "data-xlabel": [], styles: [], alpha:1 }
                );
            }

            expression_color = null
            if (SCORING_METHOD == 'gene') {
                if (expression[node]['gbc'] == 255) expression[node]['gbc'] = 254;
                if ( SELECTED_COLOR_SET == 'linear_purple' ) {
                    expression_color = COLORS[SELECTED_COLOR_SET][ expression[node]['gbac'] ];
                } else {
                    expression_color = COLORS[SELECTED_COLOR_SET][ expression[node]['gbc'] ];
                }
            } else if (SCORING_METHOD == 'tissue') {
                if (expression[node]['tbc'] == 255) expression[node]['tbc'] = 254;
                if ( SELECTED_COLOR_SET == 'linear_purple' ) {
                    expression_color = COLORS[SELECTED_COLOR_SET][ expression[node]['tbac'] ];
                } else {
                    expression_color = COLORS[SELECTED_COLOR_SET][ expression[node]['tbc'] ];
                }
            } else if (SCORING_METHOD == 'dataset') {
                if (expression[node]['dbc'] == 255) expression[node]['dbc'] = 254;
                if ( SELECTED_COLOR_SET == 'linear_purple' ) {
                    expression_color = COLORS[SELECTED_COLOR_SET][ expression[node]['dbac'] ];
                } else {
                    expression_color = COLORS[SELECTED_COLOR_SET][ expression[node]['dbc'] ];
                }
            }
            groups[group_name]['styles'].push(expression_color);
        } // end for expression

        series_data = [];
        for (var group_idx=0, group_n=group_names.length; group_idx < group_n; group_idx++) {
            group_name = group_names[group_idx]

            for (var i=0; i < max_series_size; i++) {
                if (!groups[group_name]['values']) {
                    series[i]['values'].push(null);
                    series[i]['styles'].push(null);
                    series[i]['errors'].push(null);
                    series[i]['data-xlabel'].push(null);
                    series[i]['p_values'].push(null);
                    series[i]['text'] = null;
                } else {
                    series[i]['values'].push(groups[group_name]['values'][i]);
                    series[i]['styles'].push(groups[group_name]['styles'][i]);
                    series[i]['errors'].push(groups[group_name]['errors'][i]);
                    series[i]['data-xlabel'].push(groups[group_name]['data-xlabel'][i]);
                    series[i]['p_values'].push(groups[group_name]['p_values'][i]);
                    series[i]['text'] = y_axis_labels[i]; //add label for legend
                }
            }
        }

        // Graphset vars that depend on graph_type
        var yMinValue = null;
        var yMaxValue = null;
        var x_start = '';
        var y_start = '';
        var plot_height = '';
        var plot_width = '';
        var offset_y = null;
        var font_size = null;

        var short_yaxis_toggle = false;
        var short_yaxis_unit = null;
        var short_yaxis_separator = null;

        var plotareaOptions = {};
        var legend_options = {};

        // Zingchart options for 'bar' graphs
        if (graph_type == 'bar') {

            //Set legend - hide it
            legend_options['visible'] = false;

            //Set top margin (title dependent)
            if (plot_title) {
                plotareaOptions['margin-top'] = '25px';
            } else {
                plotareaOptions['margin-top'] = '15px';
            }

            // Set bottom margin
            if (len == 1) {
                plotareaOptions['margin-bottom'] = '60px';
            } else {
                // Places group name above math button (prevents overlap)
                plotareaOptions['margin-bottom'] = '75px';
            }

            // Set left & right margins
            //    1 & 2 subplots set 'dynamic'
            //    3 subplots set manually
            if (len < 3) {
                plotareaOptions['margin-left'] = 'dynamic';
                plotareaOptions['margin-right'] = 'dynamic';
            } else {
                if (plot == 0) {
                    plotareaOptions['margin-left'] = '30px';
                    plotareaOptions['margin-right'] = '10px';
                }
                if (plot == 1) {
                    plotareaOptions['margin-left'] = '25px';
                    plotareaOptions['margin-right'] = '15px';
                }
                if (plot == 2) {
                    plotareaOptions['margin-left'] = '25px';
                    plotareaOptions['margin-right'] = '15px';
                }
            }

            // Determine Y-axis max and min values
            if (series[0]['values']) {
                var series_len=series.length;

                if (series_len > 1) {
                    // Multiple series per graph

                    //Set up Max & Min for y-axis
                    var series_val_list = [];
                    var series_err_list = [];
                    for (var series_index=0; series_index<series_len; series_index++) {
                        // Source of idea for using $.map http://stackoverflow.com/a/22962986/2900840
                        $.map(series[series_index]['values'], function(item){
                            series_val_list.push(item);
                        });
                        $.map(series[series_index]['errors'], function(item){
                            if (item) {
                              series_err_list.push(item);
                            } else {
                              series_err_list.push(0);
                            }
                        });
                    }
                    //Adjust axis max if errors are present
                    yMaxValue = (Math.max(...series_val_list) * 1.05) + Math.max(...series_err_list);

                } else {
                    // Only 1 series per graph

                    //Set up Max & Min for y-axis
                    yMaxValue = (Math.max(...series[0]['values']) * 1.05);
                    // yMinValue = Math.min(...series[0]['values']) - (Math.min(...series[0]['values']) * 0.05);

                    var series_err_list = [];
                    for (var series_index=0; series_index<series_len; series_index++) {
                        // Source of idea for using $.map http://stackoverflow.com/a/22962986/2900840
                        $.map(series[series_index]['errors'], function(item){
                            if (item) {
                              series_err_list.push(item);
                            } else {
                              series_err_list.push(0);
                            }
                        });
                    }
                    //Adjust axis max if errors are present
                    yMaxValue += Math.max(...series_err_list);
                }
            } //end if ['values']

            y_start = '0%';
            plot_height = '100%';
            offset_y = 30;

            if (len == 1) {
                // 1 Plot
                plot_width = '100%';
                x_start = '0%';

                //Set graph start position and width
                y_start = '0%';
                plot_height = '100%';
                offset_y = 30;
            } else if (len == 2) {
                // 2 Subplots
                plot_width = '50%';
                font_size = 12;

                if (plot == 0) {
                    x_start = '0%';
                } else {
                    x_start = '50%';
                }

            } else {
                // 3 Subplots
                plot_width = '33.3%';
                font_size = 10;

                if (plot == 0) {
                    x_start = '0%';
                } else if (plot == 1) {
                    x_start = '33.3%';
                } else {
                    x_start = '66.6%';
                }

            }
        } //end options for 'bar' graphs

        if (graph_type == 'line') {
            //Set graph x-axis position and width
            x_start = '0%';
            plot_width = '100%';
            offset_y = 0;

            // Set margins
            plotareaOptions['margin-left'] = '60px';
            plotareaOptions['margin-right'] = '20px';

            //Set top and bottom margin for each subplot
            if (plot_title) {
                plotareaOptions['margin-top'] = '25px';
                plotareaOptions['margin-bottom'] = '35px';
            } else {
                if (len == 1) {
                    plotareaOptions['margin-top'] = '10px';
                    plotareaOptions['margin-bottom'] = '40px';
                }
                if (len == 2) {
                    if (plot == 0) {
                        plotareaOptions['margin-top'] = '10px';
                        plotareaOptions['margin-bottom'] = '35px';
                    }
                    if (plot == 1) {
                        plotareaOptions['margin-top'] = '5px';
                        plotareaOptions['margin-bottom'] = '40px';
                    }
                }
                if (len == 3) {
                    if (plot == 0) {
                        plotareaOptions['margin-top'] = '10px';
                        plotareaOptions['margin-bottom'] = '35px';
                    }
                    if (plot == 1) {
                        plotareaOptions['margin-top'] = '5px';
                        plotareaOptions['margin-bottom'] = '40px';
                    }
                    if (plot == 2) {
                        plotareaOptions['margin-top'] = '5px';
                        plotareaOptions['margin-bottom'] = '40px';
                    }
                }
            }

            // Determine Y-axis max and min values
            if (series[0]['values']) {
                var series_len=series.length;

                if (series_len > 1) {
                    // Multiple series per graph

                    //Set up legend - only show legend for single plots containing multiple plots
                    legend_options['visible'] = true;
                    legend_options['marker'] = {'type': 'circle'};
                    legend_options['align'] = 'center';
                    legend_options['vertical-align'] = 'bottom';

                    //Tweak plot for legend
                    plotareaOptions['margin-bottom'] = '20px';

                    //Set up Max & Min for y-axis
                    var series_val_list = [];
                    var series_err_list = [];
                    for (var series_index=0; series_index<series_len; series_index++) {
                        // Source of idea for using $.map http://stackoverflow.com/a/22962986/2900840
                        $.map(series[series_index]['values'], function(item){
                            series_val_list.push(item);
                        });
                        $.map(series[series_index]['errors'], function(item){
                            if (item) {
                              series_err_list.push(item);
                            } else {
                              series_err_list.push(0);
                            }
                        });
                    }
                    yMaxValue = (Math.max(...series_val_list) * 1.05) + Math.max(...series_err_list);
                    yMinValue = Math.min(...series_val_list) - (Math.min(...series_val_list) * 0.05) - Math.min(...series_err_list);

                } else {
                    // Only 1 series per graph

                    //Set up legend - hide it
                    legend_options['visible'] = false;

                    //Set up Max & Min for y-axis
                    yMaxValue = (Math.max(...series[0]['values']) * 1.05);
                    yMinValue = Math.min(...series[0]['values']) - (Math.min(...series[0]['values']) * 0.05);

                    //Adjust axis max and min if errors are present
                    var series_err_list = [];
                    for (var series_index=0; series_index<series_len; series_index++) {
                        // Source of idea for using $.map http://stackoverflow.com/a/22962986/2900840
                        $.map(series[series_index]['errors'], function(item){
                            if (item) {
                              series_err_list.push(item);
                            } else {
                              series_err_list.push(0);
                            }
                        });
                    }
                    yMaxValue += Math.max(...series_err_list);
                    yMinValue -= Math.min(...series_err_list);
                }
            } //end if ['values']

            // Set subplot dimensions and y-start position
            if (len == 1) {
                // 1 Plot
                plot_height = '100%';
                y_start = '0%';
                y_axis_label = plot_title;
                plot_title = null;

            } else if (len == 2) {
                // 2 Subplots
                plot_height = '50%';
                font_size = 12;

                //Reverse order to match user's original intent of plot arrangement
                // (not sure why only linegraphs need this)
                series.reverse();
                y_axis_labels.reverse();
                plot_titles.reverse();


                //Set subplot start position
                if (plot == 0) {
                    y_start = '0%';
                } else {
                    y_start = '50%';
                }
            } else {
                // 3 Subplots
                plot_height = '33.3%';
                font_size = 10;

                //Reverse order to match user's original intent of plot arrangement
                // (not sure why only linegraphs need this)
                series.reverse();
                y_axis_labels.reverse();
                plot_titles.reverse();

                //Set subplot start position
                if (plot == 0) {
                    y_start = '0%';
                } else if (plot == 1) {
                    y_start = '33.3%';
                } else {
                    y_start = '66.6%';
                }
            }
        }

        if (yMaxValue > 1000){
            short_yaxis_toggle = true;
            short_yaxis_unit = "K";
            short_yaxis_separator = ",";
        }
        graphset_list.push({
            type: graph_type,
            legend: legend_options,
            title: {
              text: plot_title,
              'font-size': font_size
            },
            height: plot_height,  // height of container chart of cover
            width: plot_width,  // width of container chart will cover
            x: x_start,         // starting position of chart. start 0, 33.3% or 66.6%
            y: y_start,
            plotarea: plotareaOptions,
            plot:{
                error: {},
                valueBox:{
                  text: "%data-xlabel",
                  color: "black",
                  placement: "bottom",
                  maxChars: 7,
                  fontAngle: 90,
                  offsetY: 7
                }
            },
            scaleX:{
                labels: group_names,
                item:{
                  offsetY: offset_y
                }
            },
            scaleY:{
                label: {
                    text: y_axis_label
                },
                minValue: yMinValue,
                maxValue: yMaxValue,
                // offsetEnd: offsetyTop,
                short: short_yaxis_toggle,
                shortUnit: short_yaxis_unit,
                thousandsSeparator: short_yaxis_separator
            },
            tooltip:{
              jsRule: "CustomFn.formatTooltip()"
            },
            series: series
        });
    } //end for loop expressions_list

    CustomFn = {};

    // Add the graphset(s) to zing_config
    var zing_config = {
        graphset: graphset_list
    };
    // Create custom tooltip to show expression val and p-val
    // Source: http://stackoverflow.com/a/37038343/2900840
    CustomFn.formatTooltip = function(p) {
        var dataset = zingchart.exec(p.id, 'getdata');
        var plotIndex = p.plotindex;
        var nodeIndex = p.nodeindex;
        var series = dataset.graphset[p.graphindex].series;

        var tooltipText = "";
        for (var i=0, len=series.length; i<len; i++) {
            // add expression val
            if ( i == plotIndex ) {
                tooltipText = series[i].values[nodeIndex] + "";

                // add p-value if included
                if ( series[i].p_values[nodeIndex] != null ) {
                    tooltipText += "\np-value: " + series[i].p_values[nodeIndex].toString() + "";
                }
                // add std_dev if included (Zingchart errors are array. Check first val)
                if ( series[i].errors[nodeIndex][0] != null ) {
                    tooltipText += "\nstd dev: " + series[i].errors[nodeIndex].toString() + "";
                }
            }
        }
        // return text and css
        return {
            text: tooltipText,
            backgroundColor: "#fff",
            color: "#000",
            fontSize: "14px",
            borderColor: "black",
            borderWidth: "1px",
            borderRadius: "3px",
            padding: "3px",
        }
    }

    //Finally, Render the graph
    zingchart.render({
        id : graph_container_target_id,
        data : zing_config,
        height: 350
    });
}

function resize(){
    $("canvas").outerHeight(365 - $("canvas").offset().top - Math.abs($("canvas").outerHeight(true) - $("canvas").outerHeight()));
}

function show_search_result_info_box() {
    if( permalinked_dataset_id ) {
        // show links_out and gene_annot with zoom_on
        $('#links_out_c, #gene_details_c').addClass('search_result_c').removeClass('search_result_c_DISABLED').show('fade', {}, 400);
        $('#dataset_zoomed_c').show('fade', {}, 400);
    } else {
        // $('div.search_result_c').show('fade', {}, 400);
        $('.search_result_c_DISABLED').addClass('search_result_c').removeClass('search_result_c_DISABLED');
        // $('div.search_result_c').toggleClass('search_result_c_DISABLED');
    }
}

$('#search_results').on("click", "a", function(e) {
    e.preventDefault(); //prevent page scrolling to top
    $(this).blur(); //removes focus so active's purple coloring can show
    select_search_result(this);
});

// Warn user if no datasets in profile
$( "#search_gene_symbol").focus(function(){
    if (dataset_collection_panel.datasets.length < 1) {
        $("#search_gene_symbol").popover('show');
    } else {
        $("#search_gene_symbol").popover('hide');
    }
});
$("#search_gene_symbol").blur(function() {
    $("#search_gene_symbol").popover('hide');
});

// Popover for warning user that profile lacks datasets.
$('#search_gene_symbol').popover({
  	animation: true,
  	trigger: 'manual',
    container: 'body',
  	content: "<div class='text-center' style='width:250px;'>" +
        "<div class='alert alert-warning text-center'>" +
        "<p><span class='fa fa-exclamation'></span> <b>No datasets in current profile</b></p>" +
      	"</div>" +
        "<p>To search a gene, add a dataset to your layout profile in the <a href='./dataset_manager.html'>Dataset Manager</a>.</p>" +
        "</div>",
  	html: true,
  	placement: 'right'
});

$( "#gene_search_form" ).submit(function( event ) {
    $("#viewport_intro").hide();
    $("#viewport_main").show();

    // determine if searching for exact matches
    $("#exact_match").val( $('#exact_match_input').prop("checked") );

    $('#recent_updates_c').hide();
    $('#searching_indicator_c').show();

    var formData = $("#gene_search_form").serializeArray();
    var URL = $("#gene_search_form").attr("action");

    $('#search_results').empty();
    // show search results
    $('#search_results_c').removeClass('search_result_c_DISABLED');

    $.ajax({
        url : './cgi/search_genes.py',
        type: "POST",
        data : formData,
        dataType:"json",
        success: function(data, textStatus, jqXHR) {
        	// reset search_results
        	search_results = data;

            populate_search_result_list(data);
            $('#searching_indicator_c').hide();
            $('#intro_content_c').hide('fade', {}, 400, function() {
                // auto-select the first match.  first <a class="list-group-item"
                first_thing = $('#search_results a.list-group-item').first();
                select_search_result(first_thing);
            });

            // http://manos.malihu.gr/jquery-custom-content-scroller/
            // The author of this utility was wonderfully responsive and helpful
            if (SCROLLBAR_DRAWN == false) {
                $("#search_results_scrollbox").mCustomScrollbar({
                	theme: '3d-thick-dark',
                    scrollButtons:{ enable:true },
                    // we need to disable the keyboard scrolling so our custom indicators can work
                    keyboard:{ enable:false }
                });

                // Change height after initializing scroller.
                $("#search_results_scrollbox").css({'height': 'calc(90vh - 260px)'});

                SCROLLBAR_DRAWN = true;
            } else {
                $("#search_results_scrollbox").mCustomScrollbar("update");
            }

            return false;
        },
        error: function (jqXHR, textStatus, errorThrown) {
            $('#searching_indicator_c').hide();

            // Error occurred
      			if ( $('#search_gene_symbol').val().length < 1 ) {
                // No gene symbol entered
        				$('.alert-container').html('<div class="alert alert-danger alert-dismissible" role="alert">' +
        					'<button type="button" class="close close-alert" data-dismiss="alert" aria-label="Close"><span aria-hidden="true">&times;</span></button>' +
        					'<p class="alert-message"><strong>Oops! </strong>No gene symbol was entered. Enter a gene symbol and try again.</p></div>').show();

      			} else if ( dataset_ids_loaded.length == 0) {
                // No datasets in current layout profile
        				$('.alert-container').html('<div class="alert alert-danger alert-dismissible" role="alert">' +
        					'<button type="button" class="close close-alert" data-dismiss="alert" aria-label="Close"><span aria-hidden="true">&times;</span></button>' +
        					'<p class="alert-message"><strong>Oops! </strong>No datasets were found in the current layout profile.</p><p>To add datasets to a profile or choose a different profile, go to the <a href="./dataset_manager.html" class="alert-link">Dataset Manager</a>.</p></div>').show();

            } else {
                // Some other error occurred
              	display_error_bar(jqXHR.status + ' ' + errorThrown.name);
      		}
        }
    });

    return false;  // keeps the page from not refreshing

});

// controls to enable user scrolling of results with mouse arrow
scrolling_results = false

$('body').click(function(event) {
    if (!$(event.target).closest('#search_results_c').length) {
        scrolling_results = false
    } else {
        scrolling_results = true
    };
});

$(document).keydown(function(event) {
    // don't do anything unless we're scrolling results
    if (scrolling_results == true) {
        // this makes sure the browser doesn't scroll the window
        event.preventDefault();

        switch (event.keyCode) {
            // up key
            case 38:
            if (AT_FIRST_MATCH_RECORD == false) {
                select_search_result($(SELECTED_GENE).prev())
            }
            break;

            // down key
            case 40:
            if (AT_LAST_MATCH_RECORD == false) {
                select_search_result($(SELECTED_GENE).next())
            }
            break;
        }
    }
});

// Apply math tranformation to a dataset
$(document).on('click', '.select_math_transform', function(e) {
    e.preventDefault();

    var math_selected = $(this).data('math');
    var dataset_id = $(this).data('dataset_id');

    // if already displayed, don't apply math
    if ( math_selected === $('#math_menu_' + dataset_id).attr('data-math') ) {

    } else {
        var current_gene_sym = SELECTED_GENE.text();

        //fade dataset container until success/error
        //$('#dataset_' + dataset_id).css('opacity', '0.5');

        // get index of gene. This will help locate corresponding expression values
        var index = '';
        for (var i=0, search_len = search_results.length; i < search_len; i++) {
            var annotations = search_results[i]['annotations'];
            for (var n=0, annot_len = annotations.length; n < annot_len; n++) {
                if (annotations[n].gene_symbol === current_gene_sym) {
                    index = i;
                    break;
                }
            }
        }
        var expression = search_results[index]['expression'];

        //do math transformations
        // math results are stored beside the dataset's expression values.
        // this preserves the original raw_vals and std_dev from which all math is performed on
        for (var i=0, len=expression.length; i<len; i++) {

            //Only transform the selected dataset's values and std_dev (if present)
            if ( expression[i]['dataset_id'] == dataset_id ) {
                var raw_val = expression[i]['raw_val'];
                var std_dev = expression[i]['std_dev'];

                if ( math_selected == 'log2' ) {
                    var new_val = parseFloat(Math.log2(raw_val).toFixed(4));
                    expression[i]['new_val'] = new_val;
                    expression[i]['new_std_dev'] = null;
                } else if ( math_selected == 'log10' ) {
                    var new_val =  parseFloat(Math.log10(raw_val).toFixed(4));
                    expression[i]['new_val'] = new_val;
                    expression[i]['new_std_dev'] = null;
                } else {
                    //keep the same val
                    expression[i]['new_val'] = raw_val;
                    expression[i]['new_std_dev'] = std_dev;
                }
            } //end if dataset_id
        } //end for

        //if user is logged in
        if (session_id) {
            $.ajax({
                url: './cgi/change_dataset_user_math_preference.cgi',
                type: 'POST',
                data: { 'session_id': session_id, 'dataset_id': dataset_id, 'math_preference': math_selected },
                dataType: 'json',
                success: function(data, textStatus, jqXHR) {
                    if (data['success'] == 1) {

                      //redraw graph
                      if (dataset_id in bargraphs || dataset_id in linegraphs) {
                          // Get only dataset's expression data
                          var dataset_expression = $.grep(expression, function(a) {
                              return a.dataset_id === dataset_id;
                          });
                          draw_graphs(dataset_expression, dataset_id);
                      }

                      //With success, change the corner-badge
                      $('#math_menu_' + dataset_id).removeClass('data-raw data-log2 data-log10').addClass('data-' + math_selected)
                        .text(math_selected).attr('data-math', math_selected);
                    } else {
                        console.log(data['error']);
                    }
                },
                error: function(jqXHR, textStatus, errorThrown) {
                    $('#dataset_' + dataset_id).css('opacity', '1.0');

                    display_error_bar(jqXHR.status + ' ' + errorThrown.name);

                }
            });//end ajax
        } else {
            //redraw graph
            if (dataset_id in bargraphs || dataset_id in linegraphs) {
                // Get only dataset's expression data
                var dataset_expression = $.grep(expression, function(a) {
                    return a.dataset_id === dataset_id;
                });
                draw_graphs(dataset_expression, dataset_id);
            }

            //change the corner-badge
            $('#math_menu_' + dataset_id).removeClass('data-raw data-log2 data-log10').addClass('data-' + math_selected)
              .text(math_selected).attr('data-math', math_selected);
        } //end if/else session_id
    } // end if math is same or different than displayed
});

// Toggle plot type for a dataset
$(document).on('click', '.change_plot_type', function(e) {
    console.log("WARNING: Called deprecated function click on .change_plot_type");
    return false;

    e.preventDefault();

    var plottype_selected = $(this).data('plot-type');
    var dataset_id = $(this).data('dataset_id');

    // if already displayed, do nothing
    if ( plottype_selected === $('#plot_type_menu_' + dataset_id).attr('data-plot-type') ) {

    } else {
        var current_gene_id = SELECTED_GENE.data('id');

        // Hide old plot
        $("#dataset_" + dataset_id + "_h5ad").hide().empty();
        $('#plot_type_menu_' + dataset_id).hide();

        $.ajax({
            url: './cgi/get_h5ad_plot.cgi',
            type: 'POST',
            dataType: 'json',
            data: {'gene_id': current_gene_id, 'dataset_id': dataset_id, 'plot_format': plottype_selected},
            success: function(data, textStatus, jqXHR) {
                // Data found for gene. Show plot
                if (data['success'] == 1) {

                    $('#dataset_' + dataset_id + '_h5ad').html(data['html']);
                    $('#plot_type_menu_' + dataset_id).text(plottype_selected);

                    //Show new plot
                    $("#dataset_" + dataset_id + "_h5ad").show();
                    $('#plot_type_menu_' + dataset_id).show();
                    $("#" + dataset_id + "_dataset_status_c").hide();

                    //Resize plot
                    resizePlot( $("#dataset_" + dataset_id + "_h5ad") );
                }

                // No data for gene.
                if (data['success'] == -1) {
                    $("#" + dataset_id + "_dataset_status_c h2").text('Init.  Change me');
                }
            },
            error: function(jqXHR, textStatus, errorThrown) {
                console.log(jqXHR);
                console.log("textStatus: " + textStatus);
                console.log("errorThrown: " + errorThrown);
                display_error_bar(jqXHR.status + ' ' + errorThrown.name);
            }
        }); //end $.ajax


        //if user is logged in
        if (session_id) {
            $.ajax({
                url: './cgi/change_dataset_user_plot_preference.cgi',
                type: 'POST',
                data: { 'session_id': session_id, 'dataset_id': dataset_id, 'plot_preference': plottype_selected },
                dataType: 'json',
                success: function(data, textStatus, jqXHR) {
                    if (data['success'] == 1) {
                        console.log('new plot preference saved')
                    } else {
                        console.log(data['error']);
                    }
                },
                error: function(jqXHR, textStatus, errorThrown) {
                    display_error_bar(jqXHR.status + ' ' + errorThrown.name);
                }
            });//end ajax
        } //end if/else session_id
    } // end if plottype is same or different than displayed
});

// Gene Details collapse css changes
$(document).on('click', '#gene_details_header, #gene_collapse_btn', function() {
    if ($('#gene_collapse_btn').attr('aria-expanded') == 'false') {
        //Details is collapsed
        $('#gene_details_header').css({
            'border-bottom-left-radius': '4px',
            'border-bottom-right-radius': '4px'
        });

        //Change '-' back to '+' button
        $('#gene_collapse_btn').replaceWith('<span id="gene_collapse_btn" class="fa ' +
            'fa-minus pull-right" title="Show gene information" ' +
            'data-toggle="collapse"' +
            'data-target="#gene_details_info" aria-expanded="true" ' +
            'aria-controls="gene_details_info"></span>');
    } else {
      //Details is expanded
      $('#gene_details_header').css({
          'border-bottom-left-radius': '0px',
          'border-bottom-right-radius': '0px'
      });

      //Change '+' to '-'
      $('#gene_collapse_btn').replaceWith('<span id="gene_collapse_btn" class="fa ' +
          'fa-plus pull-right" title="Hide gene information" ' +
          'data-toggle="collapse" ' +
          'data-target="#gene_details_info" aria-expanded="false" ' +
          'aria-controls="gene_details_info"></span>');
    }
});

if (window.location.href.indexOf("manual.html") === -1) { //Without this the manual links break.
    // Allows user to change layout profile without having to go to dataset Manager
    // User not logged in can only change domain profiles
    $('#selected_profile').editable({
        mode: 'inline',
        type: 'select',
        pk: function() {
            //return the layout_id
            return $('.editable-input select option:selected').attr('value');
        },
        params: function(params) {
            //data to be sent to 'url'
            var data = {};
            data['layout_id'] = parseInt(params.pk);
            data['session_id'] =  Cookies.get('gear_session_id');
            return data;
        },
        url: function(params) {
            /*
              When a user changes the profile the following needs to happen

              - Store the new layout ID
              - Fetch and draw the new DatasetCollectionPanel
              - If a user is logged in, set a cookie with the new layout ID
             */
            dataset_collection_panel.set_layout(params.layout_id, $('.editable-input select option:selected').text(), true);
            update_datasetframes_generesults();
        },
        title: 'Click to change',
        tpl: '<select style="width:160px"></select>',
        showbuttons: false,
        // source: layouts, //This should work. Not sure why it doesn't
        source: function() {
            return layouts; //populated by load_layouts()
        },
        sourceCache: true
    });

    // Load user's gene carts
    $('#selected_gene_cart').editable({
        mode: 'inline',
        type: 'select',
        pk: function() {
            //return the gene_cart_id
            return $('.editable-input select option:selected').attr('value');
        },
        params: function(params) {
            //data to be sent to 'url'
            var data = {};
            data['gene_cart_id'] = parseInt(params.pk);
            data['session_id'] =  Cookies.get('gear_session_id');
            return data;
        },
        url: function(params) {
            var d = new $.Deferred(); //Causes editable to wait until results are returned

            //User is not logged in
            if (!params.session_id) {
                d.resolve();
            } else {
                //User is logged in
                $("#search_gene_symbol").prop("disabled", true);
                $("#selected_gene_cart_loading_c").show();

                //Get the gene cart members and populate the gene symbol search bar
                $.ajax({
                    url: './cgi/get_gene_cart_members.cgi',
                    type: 'post',
                    data: params,
                    success: function(data, newValue, oldValue) {
                        if (data['success'] == 1) {
                            // Append gene symbols to search bar
                            gene_symbols = ''

                            //format gene symbols into search string
                            $.each(data['gene_symbols'], function(i, item){
                                gene_symbols += item['label'] + ' ';
                            });

                            $('#search_gene_symbol').val(gene_symbols);

                            // determine if searching for exact matches
                            if ( $('#exact_match_input').prop("checked") == false ) {
                                $('#exact_match_input').bootstrapToggle('on');
                            }

                        } else {
                            $('#selected_gene_cart').text(oldValue);
                            $('.alert-container').html('<div class="alert alert-danger alert-dismissible" role="alert">' +
                              '<button type="button" class="close close-alert" data-dismiss="alert" aria-label="Close"><span aria-hidden="true">&times;</span></button>' +
                              '<p class="alert-message"><strong>Oops! </strong> ' + data["error"] + '</p></div>').show();
                        }
                        $("#search_gene_symbol").prop("disabled", false);
                        $("#selected_gene_cart_loading_c").hide();
                        d.resolve();
                    }
                });
            }
            return d.promise();
        },
        title: 'Click to change',
        tpl: '<select style="width:160px"></select>',
        showbuttons: false,
        source: function() {
            return gene_carts; //populated by load_gene_carts()
        },
        sourceCache: true
    });
}

// automatically reloads dataset grid and resubmits gene search
function update_datasetframes_generesults() {
    function resubmit_gene_search() {
        $('#gene_search_form').trigger('submit');
    }

    $.when( resubmit_gene_search() ).done(function(){
        // auto-select the first match.  first <a class="list-group-item"
        first_thing = $('#search_results a.list-group-item').first();
        select_search_result(first_thing);
    });
}


// MANUAL RELATED BELOW
$(document).on('click', '.manual_link', function() {
    // remove active styling from any active sidebar sections
    $('.nav-pills > li').removeClass('active');

    //get file name of section to open
    var section_to_open = $(this).attr('data-section');

    //add styling to newly active sidebar section
    //$(this).parent().addClass('active'); //didnt work if link clicked from introduction list
    $('#open_manual_' + section_to_open).parent().addClass('active');

    //hide current sections
    $('#manual_opening_intro').hide();
    $('.container_manual').hide();

    //load selected manual section
    $('#manual_content').load('./include/manual/' + section_to_open + '.html');
});

// runs when img element is loaded
function resizeMap() {
  	// resizes the area coords to match image's dynamic width and height
  	$('map').imageMapResize();
}

function showDetails(map_area) {
    // exit if outside mapped area coords
    if ( map_area.length < 0 ) {
        return
    }

    // open tooltip for currently hovered area
    //$('[data-toggle="tooltip"]').tooltip();
    $('.text_tooltip').tooltip(({
    		html: false,
    		template:  '<div class="tooltip area-ttip" role="tooltip">' +
    					'<div class="tooltip-arrow area-ttip-arrow"></div>' +
    					'<div class="tooltip-inner area-ttip-inner"></div>' +
  				  '</div>'
    		}));
    $('.html_tooltip').tooltip(({
    		html: true,
    		template: '<div class="tooltip area-ttip" role="tooltip" style="min-width: 500px; max-width: 100%;">' +
    					'<div class="tooltip-arrow area-ttip-arrow"></div>' +
    					'<div class="tooltip-inner area-ttip-inner" style="min-width: 500px; max-width: 100%;"></div>' +
  				  '</div>'
    		}));

    // get mouse coords, offset them, then apply them to tooltip position
    //http://stackoverflow.com/a/12510161/2900840
    $('[data-toggle="tooltip"]').mousemove(function(e){
        var left = e.pageX - $(this).offset().left + 30;
        var top = e.pageY - $(this).offset().top + 30;

        // tooltip for note's create button runs off screen
        if ( $(this).attr('alt') == 'note_create' ) {
            top = top - 60;
        }
        // apply to only imagemap tooltips
        $('.area-ttip').css({ top: top, left: left });
    });
}

function set_exact_match(mode) {
    if (mode == 'on') {
        $('#exact_match_input').bootstrapToggle('on');
        $("#exact_match_icon img").attr("src", "img/arrow_target_selected.png");
        $("#exact_match_icon img").attr('data-original-title', "Exact match (currently on)").tooltip('show');
    } else if (mode == 'off') {
        $('#exact_match_input').bootstrapToggle('off');
        $("#exact_match_icon img").attr("src", "img/arrow_target_unselected.png");
        $("#exact_match_icon img").attr('data-original-title', "Exact match (currently off)").tooltip('show');
    }
}
