fn1 = '/Users/khairyk/VPCZ1_kk_share/mwork/kktoolbox_local/tissue_mechanics_kktoolbox/Fly_embryo_Berkeley_stage5_last_Lmax80_with_dorsal_L56.shp3';
fn2 = '/Users/khairyk/VPCZ1_kk_share/mwork/kktoolbox_local/tissue_mechanics_kktoolbox/vfi_result_200k.shp3';

s1 = shp_surface;
s1 = s1.read_shp_surface_ascii(fn1);
s1 = set_center(s1,[0 0 0]);
m1 = get_mesh(s1);
m1 = m1.translate_to_center_of_mass();

s2 = shp_surface;
s2 = s2.read_shp_surface_ascii(fn2);
s2 = set_center(s2,[0 0 0]);
m2 = get_mesh(s2);
m2 = m2.translate_to_center_of_mass();
plot_difference(m2,m1);