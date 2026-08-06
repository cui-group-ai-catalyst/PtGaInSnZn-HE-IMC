"""Round-particle geometry (EDIT THESE NUMBERS if you know better positions).

This is the analysis region used by make_3d.py and classify_atoms.py.
If you prefer to mark the particle yourself, either:
  * edit CX / CY / R_NM below and re-run, or
  * tell me the numbers (e.g. "centre (700, 720) px, layer radii ...") and I
    will update them.
CX, CY are the particle centre in pixels of the 1315x1315 images.
R_NM[z] is the particle radius (nm) in layer z+1.
Derived from the sharpest middle layers: bright-atom centroid = (712, 723).
"""
PX_NM = 0.01138848395   # confirmed by OME PhysicalSizeX in the original RGB TIFFs
DZ_NM = 0.5
CX, CY = 712.0, 723.0                 # particle centre (px)
# Radius per layer (nm).  Crystalline body L1-11; from L12 the liquid-metal cap
# tapers smoothly (L12 a little smaller than L11, then progressively smaller).
R_NM = [2.50, 2.70, 3.05, 3.40, 3.40, 3.50, 3.75, 3.75,
        3.75, 3.60, 3.75, 3.40, 2.50, 1.80, 1.30, 0.85]
