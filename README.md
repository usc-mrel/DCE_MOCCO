# DCE_MOCCO
Direct model-based reconstruction of TK maps for accelerated DCE-MRI using a flexible MOCCO approach.

## Demo data
please download demo phantom data from: <br /> https://drive.google.com/file/d/0B4nLrDuviSiWT3ZKUmd0YjRwUEU/view?usp=sharing <br /> 
please download demo in-vivo data from: <br /> https://drive.google.com/file/d/0B4nLrDuviSiWXzJhLWFwN1c1ZG8/view?usp=sharing 

## Demo scripts
**phantom_etofts_demo.m:**
Read pre-calculated eTofts TK maps and generated k-space (ref Yannick's work), and perform MOCCO to reconstruct TK maps from under-sampled data.
Option to select different TK solver: 1. Third-party Rocketship. 2. In-house gradient solver.

**AIF_TK_patlak_demo.m:**
Read in-vivo DCE MRI data, and retrospective under-sample the data by GOCART. Perform MOCCO to jointly reconstruct both AIF and patlak TK maps from under-sampled data. 

## Functions:

**conc2Ktrans_Y.m**: Backward modeling to convert contrast concentration to TK parameter maps.<br /> 
**conc2sigD.m**: Convert contrast concentration to signal (image difference).<br /> 
**genRGA.m**: Generate randomized golden-angle radial sampling pattern.<br /> 
**model_extended_tofts_s.m**: Forward modeling from eTofts TK maps to contrast concentration.<br /> 
**Ktrans2conc.m**: Forward modeling to convert Patlak TK maps to contrast concentration.<br /> 
**sig2conc2D.m**: Convert signal (image difference) to contrast concentration.<br /> 
**multi_disp_e.m**: Utility function to visualize eTofts TK parameters.<br /> 
**CG_recon.m**: CG reconstruction of signal (image difference) from under-sampled k-space. <br /> 
**SAIF_p.m**: Generate population-averaged AIF.<br /> 
