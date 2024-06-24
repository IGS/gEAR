offset = 0

for batch_num in range(batches + 1):
    log.debug(f"Processing batch number: {batch_num} batch size: {dumper_batch_size}")
    offset = batch_num * dumper_batch_size
    log.debug(f"Writing samples between {offset} and {offset + dumper_batch_size}")

    # The files are written in a program subdirectory, so check if it exists,
    # if not create it
    output_path = output_dir + "/" + program_obj.name + "/" + grant_obj.prj_id + "--" + cohort_obj.coh_id + "--" + dul.access + "--samples-0" + str(batch_num) + ".json"

    # Attempt to make the output directory if it does not exist
    out_dir_path = Path(output_path).resolve().parent
    Path(out_dir_path).mkdir(parents=True, exist_ok=True)

    # Test that we have privileges to write the output file
    log.info(f"Writing output to: {output_path}")
    fh = open(output_path, 'w')
    fh.close()

    samples_to_write = samples_to_process[offset:offset + dumper_batch_size]

    # The samples we have loaded so far do not have all the associations, so reload the samples
    # with the associations
    samples_to_write = [sample_util.get_sample({"id": sample.id}, 
                                                assoc = ["attributes", "anatomies", "sample_assoc_sample_parent", "sample_assoc_sample_child", "sbj_ids"]) 
                                                for sample in samples_to_write]

    # Serialize the sample so it can be written out
    samples_to_write = [sample.to_serializable() for sample in samples_to_write] 

    # Serializing json
    # Writing to sample.json
    with open(output_path, "w") as outfile:
        for sample in samples_to_write:
            outfile.write(json.dumps(sample))
            outfile.write("\n")    
    log.info(f"Wrote output to: {output_path}")
