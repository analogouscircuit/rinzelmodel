# Tools for the study of synchrony enhancing neurons and their role in pitch perception

The chief reason for making this repo public is to make available the code used to produce Figure 9
in Dahlbom and Braasch (2020), shown below.  The script that runs the trials and produces the figure itself is `scripts/JASAMRStaggered.jl`.

![GitHub Logo](/figures/JASA_Figure_9.png)


The repo additionally contains much additional apparatus for experimenting with point neuron models of the synchrony-enhancing variety,
in particular those suggested by Rothman and Manis (2003) and Meng et al. (2012).  The repo also contains most of the apparatus 
necessary to reimplement the model of Huang and Rinzel (2016), though their original matlab code has since been 
[made available](https://github.com/hcc11/PitchModel). This has all been included in case anyone has interest in this material, but,
as always, this is hacked-together research code and is offered with no warranties of any kind.  Feel free to ask questions.

## References
Dahlbom, D. A. and Braasch, J. (2020). "How to pick a peak: Pitch and peak shifting in temporal models of 
pitch perception," Journal of the Acoustical Society of America **147**, 2713-2727.

Huang, C. and Rinzel, J. (2016). "A neuronal network model for pitch selectivity and representation," 
Frontiers in Computational Neuroscience **10**:57.

Meng, X., Guguet, G., and Rinzel, J. (2012). "Type III excitability, slope sensitivity, and 
coincidence detection," Discrete and Continuous Dynamical Systems Series A **32**, 2729-2757.
