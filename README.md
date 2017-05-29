# Tradict
High fidelity inference of transcriptional states using a small, machine-learned subset of the transcriptome. 

## Quick Start

See demo_tradict.m to get started on a full scale example.

## Description
Tradict is a tool for learning a subset of marker genes that statistically represent the euykaryotic transcriptome. In prospective scenarios these markers can be measured and used to reconstruct the expression of the major transcriptional programs of the cell. In our work, transcriptional programs represent the major biological processes and pathways of the cell. 

Tradict has two modes of usage. The first is **training**, in which the transcriptome is encoded into a subset of marker genes and the quantitative relationships between these genes, the transcriptional programs, and the remaining non-marker genes are learned. The second mode is **prediction**, in which given expression measurements of the markers, Tradict uses the model learned during training to reconstruct expression of transcriptional programs and the remaining non-marker genes. We refer users to the Materials and Methods section in the Supplemental Information of our paper for a more detailed description of how exactly this is done [1]. 

## Usage
Here we provide a general description of the interface met by `tradict_train` and `tradict_predict`, the training and prediction functions of Tradict. Tradict is currently implemented in MATLAB, but we are working on a R package to make Tradict more accessible and open source. In the MATLAB terminal type, `help tradict_train.m` and `help tradict_predict.m` to see exact usage details. 

`demo_tradict.m` provides a full scale example of how to run Tradict. Users can refer to this template of commands to perform their own analyses. 

### Training (encoding)
Training a Tradict model requires 4 inputs:

1. A #-training-samples x #-genes expression matrix of training transcriptomes. The units of this matrix should be 'total number of measured transcripts'.
2. A #-training-samples x 1 vector of sequencing depths, in millions of reads, for each training sample.
3. A #-genes x 1 vector of gene or transcript IDs. 
4. A #-transcriptional-programs sized set of gene sets, each of which define the constituent members of a transcriptional program. For example, for a given transcriptional program, the corresponding gene set will specify all gene/transcript IDs of the members contained within that program.

`tradict_train` outputs a `model` object, that contains various information about the transcriptional programs, marker selection, and model parameters that define how the expression of transcriptional programs and non-marker genes can be inferred from the expression measurements of the marker genes. 

### Prediction (decoding)
Using Tradict for prediction requires 3 inputs:

1. A #-samples x #-markers expression matrix where entry (i,j) is the total number of measured transcripts for marker j in sample i.
2. A #-samples x 1 vector of sequencing depths (in millions of reads). 
3. A model object obtained from `tradict_train`.

`tradict_predict` outputs a #-samples x #-transcriptional-programs matrix of expression values of the transcriptional programs defined during training. Additionally, `tradict_predict` outputs two #-samples x #-total-genes matrices of the predicted expression values, in transcripts per million (TPM) and log-TPM, of all genes in the transcriptome.

## References
1. Surojit Biswas, Konstantin Kerner, Paulo Jose Pereira Lima Teixeira, Jeffery L. Dangl, Vladimir Jojic, Philip A. Wigge. "Tradict enables accurate prediction of eukaryotic transcriptional states from 100 marker genes" (2017). Nature Communications. https://www.nature.com/articles/ncomms15309

Send questions/comments to surojitbiswas@g.harvard.edu
