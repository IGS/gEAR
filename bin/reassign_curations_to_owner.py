#!/opt/bin/python3

"""
Not sure how generally useful this will be.

During a conference Bea did multi-gene curations for a LOT of datasets, and this 
reassigns ownership of those curations to to the dataset owner.

It also sets these curations as the default multi-gene view for each, clearing any
previous defaults.

"""
import argparse
import os
import sys

lib_path = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'lib'))
sys.path.append(lib_path)
import geardb

CURATOR_USER_ID = 499
MG_PLOT_TYPES = ['heatmap', 'volcano', 'mg_violin', 'dotplot', 'quadrant']

def main():
    cnx = geardb.Connection()
    cursor = cnx.get_cursor()

    dataset_ids = list()

    layout_names = [
        'Adult',
        'Aging',
        'Central auditory system',
        'Chick',
        'Development',
        'Ear (diverse variety)',
        'Hair cell',
        'Hair cell and supporting cell transcriptomics and epigenomics, Mouse (Segil 2021)',
        'Hair cell reprogramming - transcriptomics and epigenetics (Segil)',
        'Hearing (default)',
        'Human disease',
        'Lateral wall',
        'Mouse models',
        'Mouse models 2',
        'Neuron'
    ]

    dataset_names = [] #get_dataset_names()

    layout_qry = """
          SELECT id, user_id, label
            FROM layout 
           WHERE user_id = 0
    """
    qry_args = []

    cursor.execute(layout_qry, qry_args)

    for (layout_id, user_id, label) in cursor:
        if label in layout_names:        
            layout = geardb.Layout(id=layout_id)
            layout.get_members()
            dataset_ids.extend(layout.dataset_ids())

    for name in dataset_names:
        dataset = geardb.get_dataset_by_title(title=name)

        if dataset:
            dataset_ids.append(dataset.id)

    print("Dataset ID count: {0}".format(len(dataset_ids)))

    ## For each datasets, grab the displays
     # If a display is owned by curator change it to be owned by the dataset author
     # Also clear any preferences for that dataset for the author and set the new one
    dc = geardb.DatasetCollection()
    dc.get_by_dataset_ids(dataset_ids)

    print("Got {0} datasets to process".format(len(dc.datasets)))

    for dataset in dc.datasets:
        print("Getting displays for dataset ID: {0}".format(dataset.id))
        displays = dataset.get_displays()
        print("Got {0} displays".format(len(displays)))

        if not len(displays):
            continue

        for display in displays:
            if display.user_id == CURATOR_USER_ID and dataset.owner_id != CURATOR_USER_ID and display.plot_type in MG_PLOT_TYPES:
                print("\tDisplay id {0} is of type {1}".format(display.id, display.plot_type))
                ## change display to be owned by dataset owner
                display.user_id = dataset.owner_id
                display.save()

                ## clear any existing preferences for this dataset for this user
                ## set this one as the new preference for multi-genes
                query = """
                        INSERT INTO dataset_preference (user_id, dataset_id, display_id, is_multigene)
                             VALUES (%s, %s, %s, 1)
                             ON DUPLICATE KEY UPDATE display_id = VALUES(display_id)
                """
                print("Executing update with ({0}, {1}, {2})".format(display.user_id, dataset.id, display.id))
                print("\tDataset title: {0}".format(dataset.title))
                cursor.execute(query, (display.user_id, dataset.id, display.id))
                
    cnx.commit()
    cursor.close()
    cnx.close()


def get_dataset_names():
    return [
        "Outer Hair Cell Response to PTS-inducing noise (PrestinCre;RiboTag;CBA mice), RiboTag RNA-seq",
        "Supporting Cell Response to PTS-inducing noise (Sox2Cre;RiboTag;CBA mice), RiboTag RNA-seq",
        "Spiral Ganglion Neurons (SGN) Response to PTS-inducing noise (scRNA-seq), UMAP display",
        "Spiral Ganglion Neurons (SGN) Response to PTS-inducing noise (scRNA-seq), violin plot display",
        "Lateral Wall Response to PTS-inducing noise (scRNA-seq), UMAP display",
        "Lateral Wall Response to PTS-inducing noise (scRNA-seq), violin plot display",
        "CD45+ Cochlear Immune Cells Response to PTS-inducing noise, scRNA-seq, UMAP display",
        "CD45+ Cochlear Immune Cells Response to PTS-inducing noise, scRNA-seq, violin plot display",
        "Lateral Wall Immune Cells Response to PTS-inducing noise, scRNA-seq, UMAP display",
        "Lateral Wall Immune Cells Response to PTS-inducing noise, scRNA-seq, violin plot display",
        "Supporting Cell and Hair Cells Specific Transcriptional Responses to Heat Shock in the Mouse Utricle Epithelium (Cunningham)",
        "4-8 weeks, mouse, RNA-seq, utricle, heat shock (Cunningham)",
        "P1, mouse, RNA-seq, cochlea, hair cells, Gentamicin-induced gene expression changes (Segil)",
        "6-7wks, mouse, RNA-seq, cochleae post 2h of 120dB noise exposure (Nishizaki)",
        "P60, mouse, RNA-seq, cochlea exposed to ambient (70dB), loud (94dB) and very loud (105dB) noise (Savas)",
        "Post-hatched 1-day-old, chicken, RNA-seq, basilar papillae treated with streptomycin (Nakagawa)",
        "P3-5, mouse, microarray, organ of Corti (OC) and spiral ganglion (SG) treated with kainate (KA) for 2h and 5h, 24h and 72h post culture vs Control (CTRL) (Stankovic)",
        "Mouse, RNA-seq, HEI-OC1 cells exposed to cisplatin vs untreated controls (Hati)",
        "P0, mouse, RNA-seq, epithelial vs neuronal, mesenchymal and vascular endothelial (Hertzano)",
        "8-10 weeks-old, mouse, microarray, middle ear epithelium (MEE) in vitro air-liquid interface (ALI) culture for 7 days (Bingle)",
        "RNASEQ, schwannoma cells exposed to electromagnetic field (EMF)",
        "embryonic stem cell line, microarray, induction of hair cells (Henrique)",
        "human, RNA-seq, dermal fibroblasts transdifferentiated towards the otic lineage (Schimmang)",
        "P1, P4, mouse, RNA-seq, cochlea, with and without Wnt activation (Dabdoub)",
        "P1-2, mouse, RNA-seq, LGR5+ progenitors after neomycin injury (Chai)",
        "4-6wks, mouse, scRNA-seq, utricle Atoh1-transduced supporting cells (Groves)",
        "P7, mouse, RNA-seq, cochlea, Notch1 conditional KO with and without b-catenin overexpression (Li)",
        "Hair cell induction from mESCs, 9 days in culture, RNA-seq, derived from control, Barhl1-mutants or mouse mutant for the Atoh1-bound Barhl1-enhancer (B-Ebox) (Huang)",
        "Adult, human, RNA-seq, vestibular sensory epithelia (Forge)",
        "Microarray, synaptic regeneration, cochlear hair cells and neurons (Stankovic's lab)",
        "Dying Tall and Short HCs, chicken (Heller)",
        "Dying Tall HCs, chicken - violin plot (Heller)",
        "Dying Tall HCs, chicken - trajectory colored by expression (Heller)",
        "Dying Tall HCs, chicken - trajectory colored by cell type (Heller)",
        "Dying Short HCs, chicken - violin plot (Heller)",
        "Dying Short HCs, chicken - trajectory colored by expression (Heller)",
        "Dying Short HCs, chicken - trajectory colored by cell type (Heller)",
        "Molecular anatomy of the chicken utricle (P7) - Trajectory Striola I, Violin Display,entrezID",
        "Molecular anatomy of the chicken utricle(P7)- Trajectory Striola I, Gene Expression,entrezID",
        "Molecular anatomy of the chicken utricle (P7) - Trajectory Striola I,States,entrezID",
        "Molecular anatomy of the chicken utricle (P7)- Trajectory Striola II, Violin Display,entrezID",
        "Molecular anatomy of the chicken utricle (P7) - Trajectory Striola II,Gene Expression,entrezID",
        "Molecular anatomy of the chicken utricle (P7) - Trajectory Striola II,States,entrezID",
        "Molecular anatomy of the chicken utricle(P7)-Trajectory Extrastriola, Violin Display,entrezID",
        "Molecular anatomy of the chicken utricle(P7)-Trajectory Extrastriola, Gene Expression,entrezID",
        "Molecular anatomy of the chicken utricle(P7)-Trajectory Extrastriola, States,entrezID",
        "scRNAseq, Atoh1+Ikzf2 overexpression in adult supporting cells, published UMAP (Liu, 2021)",
        "scRNAseq, Atoh1+Ikzf2 overexpression in adult supporting cells, published violin plot(Liu, 2021)",
        "scRNA-Seq of all principal neonatal cochlear cell groups,cell type,svg (Heller)",
        "scRNA-Seq of all principal neonatal cochlear cell groups,gate (Heller)",
        "scRNA-Seq of all principal neonatal cochlear cell groups,CellTrails (Heller)",
        "scRNA-Seq of all principal neonatal cochlear cell groups,cell type,color_code (Heller)",
        "scRNA-Seq of all principal neonatal cochlear cell groups,violin,version2 (Heller)",
        "P2, mouse, scRNA-seq, Organ of Corti Clustering (Waldhaus)",
        "mouse, scRNA-seq, Organ of Corti Violin Plots (Waldhaus)",
        "P2, mouse, scATAC-seq, Organ of Corti Transcription Factor Motif Accessibility (Waldhaus)",
        "Mouse, scRNAseq,postnatal mouse utricle,SWNE plot(Jan)",
        "4dpf, zebrafish, microarray, hair cell, mantle cell and skin (Hudspeth)",
        "5dpf, zebrafish, RNA-seq, hair cell (IP) and whole larvae input (IN) (Hertzano)",
        "Neuromast Regeneration,scRNAseq--All Time Points, line plot",
        "Neuromast Regeneration,scRNAseq--control",
        "Neuromast Regeneration,scRNAseq--0min",
        "Neuromast Regeneration,scRNAseq--0.5hr",
        "Neuromast Regeneration,scRNAseq--1hr",
        "Neuromast Regeneration,scRNAseq--3hr",
        "Neuromast Regeneration,scRNAseq--5hr",
        "Neuromast Regeneration,scRNAseq--10hr"
    ]
    
    
if __name__ == '__main__':
    main()
    
