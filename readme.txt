Demo script:
phantom_etofts_demo.m:
Read pre-calculated eTofts TK maps and generated k-space (ref Yannick's work), and perform MOCCO to reconstruct TK maps from under-sampled data.
Option to select different TK solver: 1. Third-party Rocketship. 2. In-house gradient solver.

AIF_TK_patlak_demo.m:
Read in-vivo DCE MRI data, and retrospective under-sample the data by GOCART. Perform MOCCO to jointly reconstruct both AIF and patlak TK maps from under-sampled data. 