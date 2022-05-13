"use strict";

class Gene {
    constructor ({id, ensembl_id, genbank_acc, organism_id, gene_symbol,
                  product, biotype} = {}) {
        this.id = id;
        this.ensembl_id = ensembl_id;
        this.genbank_acc = genbank_acc;
        this.organism_id = organism_id;
        this.gene_symbol = gene_symbol;
        this.product = product;
        this.biotype = biotype;
    }

}

class WeightedGene extends Gene {
    constructor ({...args} = {}, weights) {
        super(args);
        this.weights = weights || [];
    }

}