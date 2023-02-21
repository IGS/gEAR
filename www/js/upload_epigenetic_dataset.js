/*
 This script relies on the source having also included the
 common.js within this project
*/

// for ui elements
var filecount = 0;

// file upload return status
var files_done = {};

// keeps track of all file upload processes
var submit_counter = 0;

// only enables submit button after all uploads are done
function enable_submit() {
    if (submit_counter == 0) {
        $("#submit").removeClass("disabled");
    }

    return;
}

function disable_submit() {
    $("#submit").addClass("disabled");
}

window.onload=function() {
    // Set the curator link href for when the user gets that far
    var home_url = `./index.html`;
    var epiviz_panel_designer_url = `./epiviz_panel_designer.html`;
    $("#goHome").click(function() {
        window.location.replace(home_url);
    });
    $("#goEpivizPanel").click(function() {
        window.location.replace(epiviz_panel_designer_url);
    });

    // $("#continue_container").hide();
};

$("#upload_more").click(function(e){
    e.preventDefault();
    var addto = "#file" + filecount;
    // var addRemove = "#file" + (filecount);
    var addTo = $("#files").children().last();
    filecount = filecount + 1;
    var newFileContainer = '\
        <div class="fileitem" id="file'+ filecount +'"> \
            <div class="form-group"> \
                <button id="remove' + (filecount) + '" class="btn btn-danger remove-me">Remove</button> \
            </div> \
            <div class="form-group"> \
                <div class="custom-file"> \
                    <input id="file'+ filecount +'_file" type="file" class="custom-file-input">  \
                    <label class="custom-file-label" for="file'+ filecount +'_file">Choose BigWig/BigBed File</label> \
                </div> \
            </div> \
            <h3>OR</h3> \
            <div class="form-group"> \
              <label for="file'+ filecount +'_url">File URL (Must be publicly available and not an intranet location)</label> \
              <input id="file'+ filecount +'_url" type="text" class="form-control"> \
            </div> \
            <div class="form-group"> \
                <label for=file'+ filecount +'_title">File Title (used as display name in visualizations)</label> \
                <input id="file'+ filecount +'_title" type="text " class="form-control ">  \
                </div> \
            <div class="form-group"> \
                <label for=file'+ filecount +'_annotation">Add annotation</label> \
                <textarea class="form-control" id="file'+ filecount +'_annotation" rows="2"></textarea> \
            </div> \
            <div class="form-group"> \
                <label for="file'+ filecount + '_genome">Genome Information (mm10, hg19 etc.)</label> \
                <select class="form-control" id="file'+ filecount + '_genome"> \
                    <option value="mm10">mm10</option> \
                    <option value="hg19">hg19</option> \
                    <option value="hg38">hg38</option> \
                    <option value="marmoset">marmoset</option> \
                </select> \
            </div> \
            <div class="form-group"> \
                <label for="file'+ filecount +'_access">Choose File Access</label> \
                <select class="form-control" id="file'+ filecount +'_access"> \
                    <option value="1">Public</option> \
                    <option value="0">Private</option> \
                </select> \
            </div> \
            <div class="form-group"> \
                <div id="file'+ filecount +'_processing"> \
                    <h2 id="file'+ filecount +'_processing_status">No File Selected</h2> \
                    <div id="file' + filecount + '_progress"> \
                        <div class="file_upload_bar" style="width: 0%;"></div> \
                    </div> \
                </div> \
            </div> \
        </div> \
    '

    var newInput = $(newFileContainer);
    // var removeBtn = '';
    // var removeButton = $(removeBtn);
    addTo.after(newInput);
    // addTo.after(removeButton);
    $("#file" + filecount).attr('data-source',addTo.attr('data-source'));
    $("#count").val(filecount);

        // update this to remove the selected elemn not the last one
        $('.remove-me').click(function(e){
            e.preventDefault();
            var fileNum = this.id.charAt(this.id.length-1);
            var fileID = "#file" + fileNum;
            $(this).remove();
            $(fileID).remove();
        });

    add_file_handler(filecount);
});

/*jslint unparam: true */
/*global window, $ */
$(function () {
    'use strict';
    add_file_handler(filecount);
});

function add_file_handler (fcount) {
    var url = 'datasets_epigenetic/uploads/';
    $('#file' + fcount + "_file").fileupload({
        url: url,
        dataType: 'json',
        // maxChunkSize: 10000000, // 10 MB
        fail: function(e, data) {
            submit_counter--;
            enable_submit();

            $('#file' + fcount + '_processing_status').html("File Upload Error");
            // $('#file' + fcount + '_processing').find(".upload_spinner").hide();

            var msg = "Unable to upload file " + data.files[0]["name"] + "to the server. Please try uploading the file again.";
            alert_user(msg);
        },
        done: function (e, data) {
            // same as dataset_upload, validate bigwigs and bigbeds!
            console.log(data.result);
            console.log(data.textStatus);

            submit_counter--;
            enable_submit();

            if (data.result.files.length == 0) {
                $('#file' + fcount + '_processing_status').html("File Upload Error");
                // $('#file' + fcount + '_processing').find(".upload_spinner").hide();
                $('#file' + fcount).addClass("alert alert-danger")

                var msg = "Unable to upload file " + data.files[0]["name"] + "to the server. Please try uploading the file again.";
                alert_user(msg);
            } else {
                files_done[data.files[0]["name"]] = true;
                $('#file' + fcount + '_processing_status').html("File uploaded successfully!");
                // $('#file' + fcount + '_processing').find(".upload_spinner").hide();
                // upload_end[fcount] = true;
            }
        },
        progressall: function (e, data) {
            // later
            // $('#file' + fcount + '_processing').show();
            // upload_begin[fcount] = true;

            var progress = parseInt(data.loaded / data.total * 100, 10);
            $('#file' + fcount + '_progress .file_upload_bar').css(
                'width',
                progress + '%'
            );
        }
    });

    // $('#file' + fcount + '_processing').hide();
    // $('#file' + fcount + '_processing').find(".upload_spinner").hide();

    $('#file' + fcount + "_file").on('change',function(e){
        var id = this.id.split("_");
        $("#" + id[0]).find(".custom-file-label").html(e.target.files[0].name);

        // if the user changes the file
        $('#file' + fcount + '_processing_status').html("Uploading file");
        // $('#file' + fcount + '_processing').show();
        // $('#file' + fcount + '_processing').find(".upload_spinner").show();
        $('#file' + fcount + '_progress .bar').css(
            'width',
            0 + '%'
        );
        $('#file' + fcount).removeClass("alert alert-danger");


        submit_counter++;
        disable_submit();
    })
}

$("#submit").click(function(e) {
	// remove any warnings, if present
	$('.label-warning').remove();

    // disable the submit button and put the UID into the form
    $("#submit").text("Processing");
    $("#submit").attr("disabled", true);
    $("#upload_more").attr("disabled", true);
    $('#dataset_uid').val(dataset_uid);

    e.preventDefault();

    var i;
    var files = {};

    $("#files").children().each(function(i, elem) {
        var fid = elem.id;
        var file_input = $('#' + fid + "_file");
        var file_url = $('#' + fid + "_url");
        var file_annotation = $('#' + fid + "_annotation");
        var file_title = $('#' + fid + "_title");
        var file_access = $('#' + fid + "_access");
        var file_organism = $('#' + fid + "_genome");
        var dataset_uid = guid('long');
        var share_uid = guid('short');

        if (Object.keys(files).includes(dataset_uid)) {
            dataset_uid += i;
            share_uid += i;
        }

        files[dataset_uid] = {};
        files[dataset_uid]["dataset_uid"] = dataset_uid;
        files[dataset_uid]["share_uid"] = share_uid;

        // can be file or url but not both!
        var is_url = false;
        if (files_done[file_input.parent().find(".custom-file-label").html()]) {
            files[dataset_uid]["file_name"] = file_input.parent().find(".custom-file-label").html();
        } else {
            files[dataset_uid]["file_url"] = file_url.val();
            is_url = true;
        }

        files[dataset_uid]["file_annotation"] = file_annotation.val();
        files[dataset_uid]["file_organism"] = file_organism.val();
        files[dataset_uid]["file_title"] = file_title.val();
        files[dataset_uid]["file_access"] = file_access.val();
        files[dataset_uid]["session_id"] = CURRENT_USER.session_id;

        if (is_url || files_done[file_input.parent().find(".custom-file-label").html()]) {
            $.ajax({
                url: './cgi/load_dataset_epiviz.cgi',
                type: "POST",
                data : files[dataset_uid],
                dataType:"json",
                success: function(data, textStatus, jqXHR) {
                    if (!data.error) {
                        console.log("upload success! check db entries!")
                        $("#upload_form").hide();
                        $("#upload_controls").hide();
                        $("#continue_container").show();
                    } else {
                        var msg = "An upload error has occurred. Please try again.";
                        alert_user(msg);

                        //reset submit button area
                        $("#submit").text("Submit");
                        $("#submit").attr("disabled", false);
                    }
                },
                error: function (jqXHR, textStatus, errorThrown) {
                    console.log('jqXHR: ', jqXHR);
                    console.log('textStatus: ', textStatus);
                    console.log('errorThrown: ', errorThrown);

                    var msg = "Unable to finish dataset upload. Please try again.";
                    alert_user(msg);
                }
            }); //end ajax
        }
    });

    // don't know why cgi doesn't recognize an array of objects.
    // $.ajax({
    //     url: './cgi/load_dataset_epiviz.cgi',
    //     type: "POST",
    //     data : files,
    //     dataType:"json",
    //     success: function(data, textStatus, jqXHR) {
    //         if (data['success'] == 1) {
    //             console.log("upload success! check db entries!")
    //             $("#upload_form").hide();
    //             $("#upload_controls").hide();
    //             $("#continue_container").show();
    //         } else {
    //             var msg = "An upload error has occurred. Please try again.";
    //             alert_user(msg);
    //             //reset submit button area
    //             $("#submit").text("Submit");
    //             $("#submit").attr("disabled", false);
    //         }
    //     },
    //     error: function (jqXHR, textStatus, errorThrown) {
    //         console.log('jqXHR: ', jqXHR);
    //         console.log('textStatus: ', textStatus);
    //         console.log('errorThrown: ', errorThrown);
    //         var msg = "Unable to finish dataset upload. Please try again.";
    //         alert_user(msg);
    //     }
    // }); //end ajax
});

$("#submit_trackhub").click(function(e) {
    e.preventDefault();

    var hub_url = $("#trackhub_url").val();

    $("#submit_trackhub").text("Processing");
    $("#submit_trackhub").attr("disabled", true);

    $.ajax({
        url: '/trackhub',
        type: "GET",
        data : {
            "hub": hub_url
        },
        dataType:"json",
        success: function(data, textStatus, jqXHR) {
            console.log(data);
            if (!data.error) {

                console.log("upload success! check db entries!")
                var files = {};

                // now add these files to the db;
                data.data.forEach(function(dfile) {
                    console.log(dfile);

                    $("#trackhub_status").append("<p> Adding File " + dfile.name + ", URL: " + dfile.datasourceId + "</p>")

                    var dataset_uid = guid('long');
                    var share_uid = guid('short');

                    files[dataset_uid] = {};
                    files[dataset_uid]["dataset_uid"] = dataset_uid;
                    files[dataset_uid]["share_uid"] = share_uid;
                    files[dataset_uid]["file_url"] = dfile.datasourceId;
                    files[dataset_uid]["file_annotation"] = dfile.annotation;
                    files[dataset_uid]["file_organism"] = $("#trackhub_genome").val();
                    files[dataset_uid]["file_title"] = dfile.name;
                    files[dataset_uid]["file_access"] = "1";
                    files[dataset_uid]["session_id"] = CURRENT_USER.session_id;

                    // very redundant, refactor & add method to send load requests between this and #submit!
                    $.ajax({
                        url: './cgi/load_dataset_epiviz.cgi',
                        type: "POST",
                        data : files[dataset_uid],
                        dataType:"json",
                        success: function(data, textStatus, jqXHR) {
                            if (!data.error) {
                                console.log("upload success! check db entries!")
                                $("#upload_form").hide();
                                $("#upload_controls").hide();
                                $("#continue_container").show();

                            } else {
                                var msg = "An upload error has occurred. Please try again.";
                                alert_user(msg);

                                //reset submit button area
                                $("#submit_trackhub").text("Submit");
                                $("#submit_trackhub").attr("disabled", false);
                            }
                        },
                        error: function (jqXHR, textStatus, errorThrown) {
                            console.log('jqXHR: ', jqXHR);
                            console.log('textStatus: ', textStatus);
                            console.log('errorThrown: ', errorThrown);

                            var msg = "Unable to finish dataset upload. Please try again.";
                            alert_user(msg);
                        }
                    }); //end ajax
                });

            } else {
                var msg = "An upload error has occurred. Please try again.";
                alert_user(msg);

                //reset submit button area
                $("#submit_trackhub").text("Submit");
                $("#submit_trackhub").attr("disabled", false);
            }
        },
        error: function (jqXHR, textStatus, errorThrown) {
            console.log('jqXHR: ', jqXHR);
            console.log('textStatus: ', textStatus);
            console.log('errorThrown: ', errorThrown);

            var msg = "Unable to load Trackhub configuration. Please try again.";
            alert_user(msg);
        }
    }); //end ajax


});

/**
 * Generates a GUID string.
 * @returns {String} The generated GUID.
 * @example af8a8416-6e18-a307-bd9c-f2c947bbb3aa
 * @author Slavik Meltser (slavik@meltser.info).
 * @link http://slavik.meltser.info/?p=142
 */
function guid(uid_length) {
    function _p8(s) {
        var p = (Math.random().toString(16)+"000000000").substr(2,8);
        return s ? "-" + p.substr(0,4) + "-" + p.substr(4,4) : p ;
    }
    if (uid_length == 'long') {
        return _p8() + _p8(true) + _p8(true) + _p8();
    }
    if (uid_length == 'short') {
        return _p8();
    }
}
