

$(document).ready(function() {

    // place the plugin HTML content into the main plugin container (from the hidden div)
    $("#upload_dataset_plugin_c").html(
        $("#import_nemoarchive_datasets_html_c").html()
    );

    // then destroy the initial placement
    document.getElementById("import_nemoarchive_datasets_html_c").remove();
});
