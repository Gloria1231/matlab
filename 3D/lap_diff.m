function lap=lap_diff(phi)
global k2
lap=real(ifftn((k2.*fftn(phi))));
end