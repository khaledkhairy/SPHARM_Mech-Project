#ifndef SIMGAST_H
#define SIMGAST_H

#include <QMainWindow>
#include "ui_simgast.h"
#include <math.h>
#include <nlopt.h>


// Forward class declarations
//class vtkSmartPointer;
class vtkDataSetMapper;
class vtkPolyDataMapper;
class vtkPolyData;
class vtkActor;
class vtkRenderer;
class vtkCellArray;
class spherical_mesh;
class shp_surface;
class shp_shell;
class vtkAxesActor;
class vtkCubeAxesActor;
class vtkCamera;
class vtkScalarBarActor;
//class vtkContextView;
//class vtkChart;
//class vtkPlot;

class SimGast : public QMainWindow
{
	Q_OBJECT

public:
	// member variables
	bool opt_intersection_test_vm, opt_intersection_test_inner,opt_intersection_test_outer;
	double gamma_v, gamma_inter, gamma_T;
	nlopt_opt opt;
	std::vector<double> onlk, tnlk;	// normalization of coefficients to make optimization easier 
	std::vector<double> Energy_vec;
	int Lmax, gdim,gdim_fit_min, func_counter, opt_plot_step,version;
	shp_shell *s;		// the shell object
	spherical_mesh *smv, *gsmv;  // for visualization
	std::string filestr;// to store the filename and path for the surface and gene expression field shp3 file
	int current_view_type;		//0 for cut-plane, 1 for mid-plane morphology
	int tri_n_intersection_test;

	// for vtk widget 1
	 double _min_sf_1, _max_sf_1;
	 bool axis_on_1;
	 bool cube_axis_actor1_on;
	 bool show_active_region;
	 bool stop_optimization;
	 bool vm_actor_on;
    vtkRenderer* _ren1, *_ren2;
    vtkAxesActor* _axes_actor1;
	vtkCamera* camera1;
	vtkCubeAxesActor *_cube_axes_actor1;
	vtkPolyDataMapper* Mapper1_in, *Mapper1_out, *Mapper1_fill, *optimMapper, *vmMapper;
	vtkActor* Actor1_in, *Actor1_out, *Actor1_fill, *optimActor, *vmActor;
	vtkCellArray *vis_polys;
	vtkPolyData* vis_polydata, *vm_polydata;
	vtkScalarBarActor *scalarBar;
	//vtkContextView* XYview;
	//vtkChart* chart;
	//vtkPlot* curve;

	// general member variables
	double _bgrdr, _bgrdg, _bgrdb;
	double _scale_fac, cscale;
	int max_gL_max;


	// general methods
	void set_L_max(const int D_new);
	void enforce_symmetry();
	void enforce_symmetry_vfi();
	void import_shape_from_disc(std::string filestr);
	void set_gLmax(int D_new);
	void report(std::string str);
	void message(std::string str);
	void shell2polydata1();
	void shell2polydata1_active_region();
	void shell_cut2polydata1(const shp_shell *s, vtkPolyData* polydata_in,vtkPolyData* polydata_out,vtkPolyData* polydata_fill);
	void update_vtk();	// directs updating the correct view according to the value of current_view_type.
	void update_vtk1();	// updates the view with the cut-shell actors
	void update_vtk2(); // updates the view with the whole shape (mid-plane morphology)  actor
	void update_vtk3();// updates the view with the whole shape (mid-plane morphology) and scalar field is active region
	void update_vtk4();// updates the view with local curvature
	void update_all_sf();
	void update_vis_optim();// calles plain (no scalar field) update
	void update_vis_optim_hi_res();
	void update_comboBox_sf();
	void update_display_fields();
	double cscale_calc(shp_surface *_s);
	void generate_vm_actor();
	void generate_vm_actor_hi_res();
	void generate_XYplot(std::vector<double> y);

	// optimization variables and methods
	std::vector<bool> ixc, iyc, izc;// the boolean vectors that determine whether clks get fitted or not
									// used for enforcing symmetries during optimization
	std::vector<int> parmix;		// index into X_o (i.e. where the value of the current fitting parameter goes into X_o)
	std::vector<double> X_o;		// construct X_o prior to optimization start, update it in objective function and then use it to update srf_m

	int maxLfit;
	bool show_hi_res;
	std::vector<double> MCsig;	//sigma for each coefficient can be set separately
	int intersection_tests();
	int vtk_intersection_tests();
	void mc_start();
	void pso_start();
	void COBYLA_start();
	void SBPLX_start();
	void ISRES_start();
	void change_configuration(int &coord, int &cix, double &old_val);
	// NLopt methods
	static double myfunc(unsigned n, const double *x, double *grad, void *my_func_data);
	static double myfunc_unconstrained(unsigned n, const double *x, double *grad, void *my_func_data);
	static double myconstraint(unsigned n, const double *x, double *grad, void *data);

	//application utility methods
	QString mimeData_2_fileName( const QMimeData * mimeData );
	void dragEnterEvent( QDragEnterEvent * event );
	void dropEvent( QDropEvent * event );

	// constructor and destructor
	SimGast(QWidget *parent = 0, QFlag flags = 0);
	~SimGast();


public slots:
	virtual void on_checkBox_symmetry_stateChanged(int state);
	virtual void on_checkBox_symmetry_vfi_stateChanged(int state);
	virtual void on_checkBox_C_o_model_stateChanged(int state);
	virtual void on_checkBox_vm_stateChanged(int state);
	virtual void on_checkBox_inner_stateChanged(int state);
	virtual void on_checkBox_outer_stateChanged(int state);
	virtual void on_checkBox_HiRes_stateChanged(int state);
	virtual	void on_pushButton_load_starting_clicked();
	virtual	void on_pushButton_load_undeformed_clicked();
	virtual	void on_pushButton_reset_view_clicked();
	virtual	void on_pushButton_section_clicked();
	virtual	void on_pushButton_mid_plane_clicked();
	virtual	void on_pushButton_axis_toggle_clicked();
	virtual	void on_pushButton_rot_z_clicked();
	virtual	void on_pushButton_rot_y_clicked();
	virtual	void on_pushButton_rot_x_clicked();
	virtual void on_pushButton_set_undef_clicked();
	virtual	void on_pushButton_start_clicked();
	virtual	void on_pushButton_status_clicked();
	virtual	void on_pushButton_clear_clicked();
	virtual	void on_pushButton_active_region_view_clicked();
	virtual void on_pushButton_cube_axis_toggle_clicked();
	virtual void on_pushButton_stop_clicked();
	virtual void on_pushButton_curvature_clicked();
	virtual void on_pushButton_constraint_surface_clicked();
	virtual	void on_pushButton_start_PSO_clicked();
	virtual void on_pushButton_stop_PSO_clicked();
	virtual void on_pushButton_start_COBYLA_clicked();
	virtual void on_pushButton_stop_COBYLA_clicked();
	virtual void on_pushButton_start_SBPLX_clicked();
	virtual void on_pushButton_stop_SBPLX_clicked();
	virtual void on_pushButton_start_ISRES_clicked();
	virtual void on_pushButton_stop_ISRES_clicked();

	virtual void on_actionOpen_session_triggered();
	virtual void on_actionSave_session_triggered();
	virtual void on_actionExit_application_triggered();
	virtual void on_actionExport_obj_triggered();
	virtual void on_actionSave_shp3_triggered();

	virtual void on_lineEdit_D_editingFinished();
	virtual void on_lineEdit_active_region_curvature_editingFinished();
	virtual void on_lineEdit_active_region_curvature_2_editingFinished();
	virtual void on_lineEdit_stiff_fac_editingFinished();
	virtual void on_lineEdit_stiff_fac_2_editingFinished();
	virtual void on_lineEdit_poisson_editingFinished();
	virtual void on_lineEdit_young_editingFinished();
	virtual void on_lineEdit_VMscale_editingFinished();
	virtual void on_lineEdit_gdim_editingFinished();
	virtual void on_lineEdit_Lmax_editingFinished();
	virtual void on_lineEdit_gLmax_editingFinished();
	virtual void on_lineEdit_gamma1_editingFinished();
	virtual void on_lineEdit_gamma2_editingFinished();
	virtual void on_lineEdit_gamma3_editingFinished();
	virtual void on_lineEdit_gamma1_2_editingFinished();
	virtual void on_lineEdit_gamma2_2_editingFinished();
	virtual void on_lineEdit_gamma3_2_editingFinished();
	virtual void on_lineEdit_mp_cutoff_editingFinished();
	virtual void on_lineEdit_mp_cutoff_2_editingFinished();
	virtual void on_lineEdit_opt_plot_frequency_editingFinished();

	virtual void on_comboBox_sf1_currentIndexChanged(int val);
	virtual void on_comboBox_sf2_currentIndexChanged(int val);
	virtual void on_comboBox_sf3_currentIndexChanged(int val);
	virtual void on_comboBox_sf1_2_currentIndexChanged(int val);
	virtual void on_comboBox_sf2_2_currentIndexChanged(int val);
	virtual void on_comboBox_sf3_2_currentIndexChanged(int val);
private:
	Ui::SimGastClass ui;
};

#endif // SIMGAST_H
