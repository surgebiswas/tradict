# Tradict
High fidelity reconstruction of the transcriptome using a small, machine-learned subset of it. 

## Quick Start

See demo_tradict.m to get started on a full scale example.

## Description

Tradict is a tool compressing the euykaryotic transcriptome into a machine-learned subset of marker genes. In prospective scenarios these markers can be measured and used to reconstruct the expression of the major transcriptional programs of the cell as well as the expression of all the remaining non-marker genes. In our work, transcriptional programs represent the major biological processes and pathways of the cell. 

Tradict has two modes of usage. The first is **training**, in which the transcriptome is encoded into a subset of marker genes and the quantitative relationships between these genes, the transcriptional programs, and the remaining non-marker genes are learned. The second mode is **prediction**, in which given expression measurements of the markers, Tradict uses the model learned during training to reconstruct expression of transcriptional programs and the remaining non-marker genes. We refer users to the Materials and Methods section in the Supplemental Information of our paper for a more detailed description of how exactly this is done [1]. 

## Usage



## References
1. Surojit Biswas, Konstantin Kerner, Paulo Jose Pereira Lima Teixeira, Jeffery L. Dangl, Vladimir Jojic, Philip A. Wigge. "Tradict enables high fidelity reconstruction of the euykaryotic transcriptome from 100 marker genes." Submitted (2016). BioRxiv: http://biorxiv.org/content/early/2016/06/21/060111.

Send questions/comments to surojitbiswas@g.harvard.edu
