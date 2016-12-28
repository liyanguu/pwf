/* MPS name: MINQ */

/* objective function: */
min: + 1 XFOU + 1 XTHR + 1 XTWO + 1 XONE ;

/* constraints */
V4 : + 0.351 XFOU + 0.578 XTHR + 0.674 XTWO + 0.785 XONE >= 0.112;
V3 : + 0.578 XFOU + 0.578 XTHR + 0.674 XTWO + 0.785 XONE >= 0.1;
V2 : + 0.674 XFOU + 0.674 XTHR + 0.674 XTWO + 0.785 XONE >= 0.095;
V1 : + 0.785 XFOU + 0.785 XTHR + 0.785 XTWO + 0.785 XONE >= 0.102;
