1. add epiviz as a resource 
2. add epiviz display type
3. add to mysql database - dataset table
4. add a dataset_display type - dataset_display type

insert into dataset_display (dataset_id, user_id, label, plot_type, plotly_config) values ("31238h1j-0000-11e7-9560-3ca9f41a0df0", 191, "Epiviz", "epiviz", '{"genes_track": {"id": "mm10"}}');

update dataset_display set plotly_config = '{"genes_track": {"id": "mm10"}}' where dataset_id = "31238h1j-0000-11e7-9560-3ca9f41a0df0"

5. add to dataset_preference

insert into dataset_preference (user_id, dataset_id, display_id) values (191, "31238h1j-0000-11e7-9560-3ca9f41a0df0", 319)

6. get genome location from gene symbol

- organism is stored as part of dataset

7. add end point for access to data /epiviz -> api.py & epiviz_data under resource

8. js/classes -> epivizdisplay - modify as needed


## Ronna/Yang dataset for gEAR

1. add to dataset - 

INSERT INTO dataset (id,owner_id,title,organism_id,pubmed_id,geo_id,is_public,ldesc,date_added,dtype,schematic_image,share_id,math_default,marked_for_removal,load_status,has_h5ad,platform_id,instrument_model,library_selection,library_source,library_strategy,contact_email,contact_institute,contact_name,annotation_source,plot_default,annotation_release) VALUES 
('8d9d963e-cb33-451b-a947-fa3a2749d57a',191,'Epiviz track viewer',1,'26328750',NULL,0,'Ronna/Yang workspace for gEAR paper',NULL,'epiviz','img/schematic_epiviz-chipseq.png','1','raw',0,'completed',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);

2. dataset_display

insert into dataset_display (dataset_id, user_id, label, plot_type, plotly_config) values ("8d9d963e-cb33-451b-a947-fa3a2749d57a", 191, "Epiviz", "epiviz", '{"genome":"mm10","fullviewer":"http://epiviz.umgear.org/yang","dataserver":"https://epiviz.umgear.org:8820","files":[],"tracks":{"epiviz-genes-track":{"id":["mm10"]},"epiviz-stacked-blocks-track":{"id":["s1peaks","s3peaks","s2peaks","s4peaks"]},"epiviz-multistacked-line-track":{"id":["s1","s3","s2","s4"],"settings":{"title":"ATAC Signal","marginTop":25,"marginBottom":23,"marginLeft":20,"marginRight":10,"measurementGroupsAggregator":"mean-stdev","step":1,"showPoints":false,"showLines":true,"showFill":true,"showErrorBars":true,"pointRadius":1,"lineThickness":2,"yMin":0,"yMax":5,"interpolation":"basis","abLine":"default"}}}}');

3. dataset_preference
 ----  change display id here from 1.
insert into dataset_preference (user_id, dataset_id, display_id) values (191, "8d9d963e-cb33-451b-a947-fa3a2749d57a", 320)
