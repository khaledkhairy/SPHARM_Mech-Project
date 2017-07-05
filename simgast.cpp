// includes for NLopt

//#include <omp.h>

// includes for QT
#include <qapplication.h>
#include <qfiledialog.h>
#include "qmenubar.h"
#include <QApplication>
#include <QFont>
#include <QPushButton>
#include <QWidget>
#include <QSlider>
#include <QtGui>
#include <QDebug>
#include <QDateTime>
#include <QCheckBox>
#include <QProgressDialog>
#include <QProgressBar>
#include <QMessageBox>

// includes for vtk
///////////////////////////////////////////////////
//// to fix problem with vtkRenderinContextOpenGL2
//// which is needed for the XY plot
#include <vtkAutoInit.h>
VTK_MODULE_INIT(vtkRenderingContextOpenGL2);
VTK_MODULE_INIT(vtkRenderingOpenGL2)
VTK_MODULE_INIT(vtkInteractionStyle)
///////////////////////////////////////////////////
#include "vtkSmartPointer.h"
#include "vtkFloatArray.h"
#include "vtkMath.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkSmartPointer.h"
#include "vtkContextView.h"
#include "vtkViewsContext2DModule.h"
#include "vtkContextScene.h"
#include "vtkChartXY.h"
#include "vtkPlot.h"
#include "vtkPlotLine.h"
#include "vtkTable.h"
#include "vtkTimerLog.h"
#include "vtkQtTableView.h"
#include "vtkPen.h"
#include "vtkAxis.h"
#include <vtkActor.h>
#include <vtkDepthSortPolyData.h>
#include <vtkCubeAxesActor.h>
#include <vtkScalarBarActor.h>
#include <vtkAxesActor.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include "vtkCylinderSource.h"
#include <vtkDataSetMapper.h>
#include "vtkCellArray.h"
#include "vtkFloatArray.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkLookupTable.h"
#include "vtkPolyDataMapper.h"
#include "vtkCurvatures.h"
#include "vtkCleanPolyData.h"
#include "vtkCamera.h"
#include "vtkOBJExporter.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkPLYWriter.h"
#include "vtkSTLWriter.h"
#include "vtkExtractEdges.h"
#include "vtkProperty.h"
#include <vtkTextProperty.h>
#include <vtkIntersectionPolyDataFilter.h>
#include "QVTKWidget.h"

#include <time.h>
#include <ctime>
#include <fstream>
#include <istream>
#include "simgast.h"
#include "shape_tools.h"


SimGast::SimGast(QWidget *parent, QFlag flags)
	: QMainWindow(parent, flags)
{
	
	QProgressDialog progress("Starting SPHARM-MECH components...", "Abort", 0, 100, 0);
    progress.setWindowModality(Qt::WindowModal);

	ui.setupUi(this);
	std::srand(std::time(0));
	std::ostringstream os;

    
	this->version = 2000;
	this->Lmax = 36;	//36
	this->max_gL_max = 36; //70
	this->gdim = 60;		//60
	this->gdim_fit_min = 30;//60
	int calc_tri_n = 3; //3
	int vis_tri_n = 5;	//5
	tri_n_intersection_test = 3;
	
	
	//omp
	//omp_set_num_threads(4);
	//Eigen::setNbThreads(4);

    sleep(1);
	progress.setValue(10);
	progress.setLabelText("Initializing the shell object: 1 of 3 (please wait)....                     ");QApplication::processEvents();
    progress.show(); progress.raise();progress.activateWindow();QApplication::processEvents();
    //os<<"Initializing shell object: 1 of 3 (please wait)...."<<std::endl;report(os.str());os.str("");
    
	/// initialize the shell object
	this->s = new shp_shell(Lmax,gdim,calc_tri_n, max_gL_max);
	//s = new shp_shell(6,30,6);
	this->s->srf_m->update();
	this->s->srf_m->update_tri();
	
    sleep(1);
	progress.setValue(20);
	progress.setLabelText("Initializing the shell object: 2 of 3 (please wait)....                     ");QApplication::processEvents();
    progress.show(); progress.raise();progress.activateWindow();QApplication::processEvents();
    //os<<"Initializing shell object: 2 of 3 (please wait)...."<<std::endl;report(os.str());os.str("");
	this->s->srf_u->update();
	this->s->srf_u->update_tri();

    sleep(1);
	progress.setValue(30);
	progress.setLabelText("Initializing the shell object: 3 of 3 (please wait)....                     ");QApplication::processEvents();
    progress.show(); progress.raise();progress.activateWindow();QApplication::processEvents();
	//os<<"Initializing shell object: 3 of 3 (please wait)...."<<std::endl;report(os.str());os.str("");
    this->s->srf_g->update();
	this->s->srf_g->update_tri();

    sleep(1);
	progress.setValue(40);
	progress.setLabelText("Initializing the surface visualization components: 1 of 2 (please wait) ....");QApplication::processEvents();
    progress.show(); progress.raise();progress.activateWindow();QApplication::processEvents();
	//os<<"Initializing surface visualization components: 1 of 2 (please wait) ...."<<std::endl;report(os.str());os.str("");
    
    /// initialize the visualization mesh objects
	smv = new spherical_mesh(this->s->srf_m->b->L_max, vis_tri_n);
	
    sleep(1);
    progress.setValue(80);
	progress.setLabelText("Initializing the surface visualization components: 2 of 2 (please wait) ....");
    progress.show(); progress.raise();progress.activateWindow();QApplication::processEvents();
    //os<<"Initializing surface visualization components: 2 of 2 (please wait) ...."<<std::endl;report(os.str());os.str("");
    gsmv= new spherical_mesh(this->s->get_gL_max(), vis_tri_n);
	s->sf1_vis.resize(gsmv->n_points);
	s->sf2_vis.resize(gsmv->n_points);
	s->sf3_vis.resize(gsmv->n_points);
	s->sf1_2_vis.resize(gsmv->n_points);
	s->sf2_2_vis.resize(gsmv->n_points);
	s->sf3_2_vis.resize(gsmv->n_points);

    sleep(1);
	progress.setValue(90);
	progress.setLabelText("Initializing optimization components....                                 ");QApplication::processEvents();
    progress.show(); progress.raise();progress.activateWindow();QApplication::processEvents();
	//os<<"Initializing optimization components...."<<std::endl;report(os.str());os.str("");
    //initialize optimization
	int nc = this->s->srf_m->xc.rows();
	ixc.resize(nc);std::fill(ixc.begin(), ixc.end(), true);
	iyc.resize(nc);std::fill(iyc.begin(), iyc.end(), true);
	izc.resize(nc);std::fill(izc.begin(), izc.end(), true);

	this->maxLfit = this->s->srf_m->b->L_max;
	opt_intersection_test_vm = true;
	opt_intersection_test_inner = true;
	opt_intersection_test_outer = true;
	show_hi_res = true;
	opt_plot_step = 20;
	
	gamma_v = 1e10;
	gamma_inter = 1e15;
	gamma_T = 1e3;
	
	func_counter = 0;
	/// initialize general settings
	
	_scale_fac = 3;
	filestr = "";
	show_active_region = 0;
	stop_optimization = false;
	vm_actor_on = false;

    sleep(1);
	progress.setValue(92);
	progress.setLabelText("Configuring visualization environment....                                    ");QApplication::processEvents();
    progress.show(); progress.raise();progress.activateWindow();QApplication::processEvents();
	// setting up the look up table
	vtkLookupTable *lut= vtkLookupTable::New();
		lut->SetNumberOfTableValues(256);
		lut->SetHueRange(0.667,0.0);
		lut->SetSaturationRange(1,1);
		lut->SetValueRange(1,1);
	_min_sf_1 = -0.5;
	_max_sf_1 = 1.3;
	// QT/VTK ---- 
	_ren1 = vtkRenderer::New();
	ui.vtkWidget1->GetRenderWindow()->AddRenderer(_ren1);
	_ren2 = vtkRenderer::New();
	ui.vtkWidget2->GetRenderWindow()->AddRenderer(_ren2);

	current_view_type = 1;
	axis_on_1 = false;
	cube_axis_actor1_on = false;
	_cube_axes_actor1 = vtkCubeAxesActor::New();

	_axes_actor1	= vtkAxesActor::New();
	_axes_actor1->SetShaftTypeToLine();
	_axes_actor1->SetConeRadius(0.1);
	_axes_actor1->SetXAxisLabelText("x");
	_axes_actor1->SetYAxisLabelText("y");
	_axes_actor1->SetZAxisLabelText("z");
	_axes_actor1->SetTotalLength(2.5, 2.5, 2.5);
	
	s->srf_m->scale(100);
	cscale = 1.0/cscale_calc(this->s->srf_m);
	_axes_actor1->SetTotalLength(_scale_fac/cscale, _scale_fac/cscale, _scale_fac/cscale);
	s->srf_m->needsUpdating = YES;

	Mapper1_in		= vtkPolyDataMapper::New();	
	Mapper1_out		= vtkPolyDataMapper::New();	
	Mapper1_fill	= vtkPolyDataMapper::New();	
	optimMapper		= vtkPolyDataMapper::New();	
	vmMapper		= vtkPolyDataMapper::New();	
	

	Actor1_in		= vtkActor::New();
	Actor1_out		= vtkActor::New();
	Actor1_fill		= vtkActor::New();
	optimActor		= vtkActor::New();
	vmActor			= vtkActor::New();

	_bgrdr = 0.5;
	_bgrdg = 0.5;
	_bgrdb = 0.5;
	_ren1->SetBackground(_bgrdr,_bgrdg,_bgrdb);
	_ren2->SetBackground(1.0, 1.0, 1.0);
	//////////////////// Set up the XY plot
	XYview = vtkContextView::New();
	XYview->GetRenderer()->SetBackground(1.0, 1.0, 1.0);
	XYview->SetInteractor(ui.vtkWidget2->GetInteractor());
	ui.vtkWidget2->SetRenderWindow(XYview->GetRenderWindow());
	//////////////////////////////////

	vtkTextProperty* tprop = vtkTextProperty::New();
    tprop->SetColor(0.02,0.06,0.62);
    tprop->SetFontFamilyToArial();
	tprop->SetFontSize(32);

	scalarBar = vtkScalarBarActor::New();
	scalarBar->SetNumberOfLabels(4);
	scalarBar->SetTitleTextProperty(tprop);
	scalarBar->SetTitle("Curvature");
	scalarBar->SetMaximumWidthInPixels(80);
	scalarBar->SetLabelTextProperty(tprop);
	//_ren1->AddActor(scalarBar);

    sleep(1);
	progress.setValue(95);
	progress.setLabelText("Inializing fast processing capability....                                    ");
    progress.show(); progress.raise();progress.activateWindow();QApplication::processEvents();
	//// prepare vtk for fast plotting during optimization by predefining the polygon structure stored in this->smv
	vis_polys = vtkCellArray::New();
	for (int i=0; i<this->smv->n_faces; i++)
	{
		vis_polys->InsertNextCell(3);
		vis_polys->InsertCellPoint(this->smv->f0[i]-1);
		vis_polys->InsertCellPoint(this->smv->f1[i]-1);
		vis_polys->InsertCellPoint(this->smv->f2[i]-1);
	}
	vis_polydata	= vtkPolyData::New();
	vis_polydata->SetPolys(vis_polys);//polys->Delete();
	/////
	vm_polydata	= vtkPolyData::New();
	vm_polydata->SetPolys(vis_polys);//polys->Delete();
	/////
	progress.setValue(100);
	progress.setLabelText("Launching application!                                                       ");
    progress.show(); progress.raise();progress.activateWindow();QApplication::processEvents();
    sleep(1);
    this->update_vtk();
    
	/////
	//os<<*s;
	//os<<"Number of available threads: "<<Eigen::nbThreads()<<std::endl;//omp_get_num_threads()
    
    os<<"Spherical harmonics tissue mechanics modeling of morphogenesis (high resolution version)"<<std::endl<<"Created by Khaled Khairy and Philipp Keller."<<std::endl;
    os<<"Copyright 2017 Howard Hughes Medical Institute"<<std::endl<<"-------------------------------------------------------------"<<std::endl<<std::endl;
    report(os.str());os.str("");
}

SimGast::~SimGast()
{

}

// utility functions
void SimGast::set_L_max(const int D_new)
{
	
	if (D_new>0 & D_new!= this->Lmax)
	{
		std::ostringstream os;os<<"Setting new L_max value "<<D_new<<" with gdim: "<<this->gdim<<" please wait ..... ";report(os.str());os.str("");
		this->Lmax = D_new;
		this->s->set_L_max(D_new);
		//this->s->set_undeformed_shape();
		s->srf_m->needsUpdating = YES;
		s->srf_u->needsUpdating = YES;
		s->srf_m->update();
		s->srf_u->update();
		s->fix_normals_sign_known =0;
		s->fix_normals_sign = 1;
		s->update();
		spherical_mesh* tmp_smv = this->smv;
		spherical_mesh* tmp_gsmv = this->gsmv;
		smv = new spherical_mesh(this->s->srf_m->b->L_max, tmp_smv->tri_n);
		gsmv= new spherical_mesh(this->s->srf_g->b->L_max, tmp_gsmv->tri_n);
		delete tmp_smv;
		delete tmp_gsmv;
		
		s->sf1_vis.resize(gsmv->n_points);
		s->sf2_vis.resize(gsmv->n_points);
		s->sf3_vis.resize(gsmv->n_points);
		generate_vm_actor();
		cscale = 1.0/cscale_calc(this->s->srf_m);
		_axes_actor1->SetTotalLength(_scale_fac/cscale, _scale_fac/cscale, _scale_fac/cscale);
		
		if(this->ui.checkBox_symmetry->isChecked()){this->enforce_symmetry();}
		if(this->ui.checkBox_symmetry_vfi->isChecked()){this->enforce_symmetry_vfi();}

		os<<"Done!\n";report(os.str());os.str("");
		update_vtk();
	}
}
void SimGast::set_gLmax(int D_new)
{
	max_gL_max = D_new;
	s->set_gL_max(D_new);
	spherical_mesh *tmp = gsmv;
	gsmv= new spherical_mesh(D_new, tmp->tri_n);
	delete tmp;

	// update the sf_vis vectors
	int val;
	val = this->ui.comboBox_sf1->currentIndex();
	this->s->sf1_vis.clear();s->sf1_vis.resize(gsmv->n_points);
	if(val>3)
	{
		gsmv->sfGen(s->srf_g->sfmx.col(val-4));	// generates the scalar field and stores it in gsmv->sf
		this->s->sfindx[0] = val-4;
		for(int i = 0;i<gsmv->n_points;i++){s->sf1_vis[i] = gsmv->sf[i];}	// copy contents of scalar field into scalar field visualization vector
	}
	
	val = this->ui.comboBox_sf2->currentIndex();
	this->s->sf2_vis.clear();
	s->sf2_vis.resize(gsmv->n_points);
	if(val>3)
	{
		gsmv->sfGen(s->srf_g->sfmx.col(val-4));	// generates the surface and stores it in gsmv->sf
		this->s->sfindx[1] = val-4;
		for(int i = 0;i<gsmv->n_points;i++){s->sf2_vis[i] = gsmv->sf[i];}// copy contents of scalar field into scalar field visualization vector
	}

	val = this->ui.comboBox_sf3->currentIndex();
	this->s->sf3_vis.clear();
	s->sf3_vis.resize(gsmv->n_points);
	if(val>3)
	{
		gsmv->sfGen(s->srf_g->sfmx.col(val-4));	// generates the surface and stores it in gsmv->sf
		this->s->sfindx[1] = val-4;
		for(int i = 0;i<gsmv->n_points;i++){s->sf3_vis[i] = gsmv->sf[i];}// copy contents of scalar field into scalar field visualization vector
	}
}
void SimGast::enforce_symmetry()
{
	this->s->srf_m->enforce_mirror_yz(this->ixc, this->iyc, this->izc);
	ixc[0] = false;
	iyc[0] = false;
	izc[0] = false;
	this->s->srf_u->enforce_mirror_yz();
	this->update_vtk();
}
void SimGast::enforce_symmetry_vfi()
{
	/*
	this->s->srf_m->enforce_vfi_pruning(this->ixc, this->iyc, this->izc);
	ixc[0] = false;
	iyc[0] = false;
	izc[0] = false;
	this->s->srf_u->enforce_vfi_pruning();
	this->update_vtk();
	*/
}
double SimGast::cscale_calc(shp_surface *_s)
{
	// find maximum value among  the L = 1 coefficients
	double mc = 0.0; //mc = max.coeff.
	for(int i=1;i<4;i++){if(std::abs(_s->xc(i))>mc){mc = std::abs(_s->xc(i));}}
	for(int i=1;i<4;i++){if(std::abs(_s->yc(i))>mc){mc = std::abs(_s->yc(i));}}
	for(int i=1;i<4;i++){if(std::abs(_s->zc(i))>mc){mc = std::abs(_s->zc(i));}}
	return abs(mc);	// dividing by the largest value
}

void SimGast::report(std::string str)
{
	//QDateTime t=QDateTime::currentDateTime();
	//this->ui.textBrowser_report->append(t.toString("hh:mm:ss"));
	QString qstr = (const char*)str.c_str();
	this->ui.textBrowser_report->append(qstr);
	//this->ui.textBrowser_report->update();
	this->ui.textBrowser_report->repaint();
	//this->ui.textBrowser_report->verticalScrollBar()->setSliderPosition(this->ui.textBrowser_report->verticalScrollBar()->maximum());
}
void SimGast::message(std::string str)
{
	QString qstr = (const char*)str.c_str();
	QMessageBox msgBox;
	msgBox.setText(qstr);
	msgBox.setInformativeText("Press OK to continue.");
	msgBox.setStandardButtons(QMessageBox::Ok);
	msgBox.setDefaultButton(QMessageBox::Ok);
	int ret = msgBox.exec();
	
}

void SimGast::import_shape_from_disc(std::string filestr)
{
	QProgressDialog progress("Starting SimGast components...", "Abort", 0, 100, 0);
    progress.setWindowModality(Qt::WindowModal);QApplication::processEvents();
	std::ostringstream os;
		os<<"Importing shp3 data from disk.\nFile: "<<filestr.c_str()<<std::endl<<" please wait ..... ";report(os.str());os.str("");
		this->ui.comboBox_sf1->setCurrentIndex(0);
		this->ui.comboBox_sf2->setCurrentIndex(0);
		this->ui.comboBox_sf3->setCurrentIndex(0);
		//os<<"Generating scalar field basis ... ";report(os.str());os.str("");
		//this->set_gLmax(this->max_gL_max);
	progress.setValue(20);
	progress.setLabelText("Configuring gene expression pattern .....                               ");//QApplication::processEvents();
		os<<"\nConfiguring gene expression pattern ..... ";report(os.str());os.str("");
		this->s->srf_g->read_trunc(this->filestr.c_str());	// read the surface with gene expression data up to gL_max
		this->s->srf_g->center_to_zero();
		s->srf_g->update();

	progress.setValue(40);
	progress.setLabelText("Configuring undeformed geometry .....                                   ");
    //QApplication::processEvents();
		//this->max_gL_max = s->srf_g->b->L_max;
		os<<"\nConfiguring undeformed geometry ..... ";report(os.str());os.str("");
		this->s->srf_u->read_trunc(this->filestr.c_str());	// read the surface with gene expression data up to gL_max
		this->s->srf_u->center_to_zero();
		s->srf_u->needsUpdating = YES;
		s->srf_u->update();
	progress.setValue(60);
	progress.setLabelText("Configuring working geometry .....                                      ");
    //QApplication::processEvents();
		os<<"\nConfiguring working geometry ..... ";report(os.str());os.str("");
		this->s->srf_m->read_trunc(this->filestr.c_str());		// read_trunc causes truncation of series or filling with zeros to keep L_max the same

		this->s->srf_m->center_to_zero();
		s->srf_m->needsUpdating = YES;
		s->srf_m->update();
	progress.setValue(80);
	progress.setLabelText("Configuring shell material .....                                        ");
    //QApplication::processEvents();
		os<<"\nConfiguring shell object ..... ";report(os.str());os.str("");
		this->s->set_undeformed_shape();		// use srf_m as the undeformed shape srf_u
		s->fix_normals_sign_known =0;
		s->fix_normals_sign = 1;
		s->update();
		this->update_comboBox_sf();
		this->s->sfindx[0] = -1;
		this->s->sfindx[1] = -1;
		this->s->sfindx[2] = -1;
		this->s->sfindx[3] = -1;
		this->s->sfindx[4] = -1;
		this->s->sfindx[5] = -1;
	progress.setValue(90);
	progress.setLabelText("Configuring display parameters .....                                    ");
    //QApplication::processEvents();
		os<<"\nConfiguring display parameters ..... ";report(os.str());os.str("");
		generate_vm_actor();
		cscale = 1.0/cscale_calc(this->s->srf_m);
		_axes_actor1->SetTotalLength(_scale_fac/cscale, _scale_fac/cscale, _scale_fac/cscale);
		this->setWindowTitle(filestr.c_str());
		update_vtk();
		os<<"\nDone!\n";report(os.str());os.str("");
	progress.setValue(100);
		/////_ren1->GetRenderWindow()->Render();
		//this->update();
}

// material properties
void SimGast::on_lineEdit_young_editingFinished()
{
	const QString text  = this->ui.lineEdit_young->text();
	double D_new = text.toDouble();
	s->set_Young(D_new); 
	//std::ostringstream os;os<<"Setting tissue Young's modulus (in Pascal) to: "<<D_new;report(os.str());
	this->update_vtk();
}
void SimGast::on_lineEdit_poisson_editingFinished()
{
	const QString text  = this->ui.lineEdit_poisson->text();
	double D_new = text.toDouble();
	s->set_Poisson(D_new); 
	//std::ostringstream os;os<<"Setting tissue Poisson ratio to: "<<D_new;report(os.str());
	this->update_vtk();
}
void SimGast::on_lineEdit_D_editingFinished()
{
	const QString text  = this->ui.lineEdit_D->text();
	double D_new = text.toDouble();
	s->set_D(D_new);
	//std::ostringstream os;os<<"Setting new average shell thickness to "<<D_new;report(os.str());
	this->generate_vm_actor();
	this->update_vtk();
}



// scalar fields
void SimGast::on_lineEdit_stiff_fac_editingFinished()
{
	const QString text  = this->ui.lineEdit_stiff_fac->text();
	double D_new = text.toDouble();
	s->set_mp_fac(D_new); 
	//std::ostringstream os;os<<"Setting active region stiffness factor to "<<D_new;report(os.str());
	//this->update_vtk();
}
void SimGast::on_lineEdit_active_region_curvature_editingFinished()
{
	const QString text  = this->ui.lineEdit_active_region_curvature->text();
	double D_new = text.toDouble();
	s->set_mp_co(D_new); 
	//std::ostringstream os;os<<"Setting active region preferred curvature to "<<D_new;report(os.str());
	//this->update_vtk();
}


void SimGast::on_lineEdit_mp_cutoff_editingFinished()
{
	const QString text  = this->ui.lineEdit_mp_cutoff->text();
	double D_new = text.toDouble();
	s->set_mp_cutoff(D_new); 
	//std::ostringstream os;os<<"Setting cutoff for scalar field to: "<<D_new;report(os.str());
	//if(this->ui.comboBox_sf1->currentIndex()>0 & this->ui.comboBox_sf2->currentIndex()>0 & this->ui.comboBox_sf3->currentIndex()>0)
	//{
		//on_pushButton_active_region_view_clicked();
		this->s->update_scalar_fields();
		std::ostringstream os;
		os<<"\n"<<"Fold 1: \n\tArea fraction vis"<<s->active_area_fraction_vis<<"\n\tArea fraction GQ:"<<s->active_area_fraction<<"\n\tMP count vis:"<<s->MP1_count<<"\n\tTotal vis:"<<s->Total_count;
	   //os<<"\n"<<"Active region area fraction (Fold 2):"<<s->active_area_fraction_vis_2<<"\t("<<s->active_area_fraction_2<<")\n";
		report(os.str());
		this->update_vtk();

	//}
	//this->update_vtk();
}
void SimGast::on_lineEdit_gamma1_editingFinished()
{
	const QString text  = this->ui.lineEdit_gamma1->text();
	double D_new = text.toDouble();
	s->set_psnl(D_new); 
	//std::ostringstream os;os<<"Setting exponent for scalar field to: "<<D_new;report(os.str());
	//if(this->ui.comboBox_sf1->currentIndex()>0 & this->ui.comboBox_sf2->currentIndex()>0 & this->ui.comboBox_sf3->currentIndex()>0)
	//{
	this->s->update_scalar_fields();
	std::ostringstream os;
	os<<"\n"<<"Active region area fraction (Fold 1):"<<s->active_area_fraction_vis<<"\t("<<s->active_area_fraction<<")\n";
	os<<"\n"<<"Active region area fraction (Fold 2):"<<s->active_area_fraction_vis_2<<"\t("<<s->active_area_fraction_2<<")\n";
	report(os.str());		
	this->update_vtk();
	//}
	
}
void SimGast::on_lineEdit_gamma2_editingFinished()
{
	const QString text  = this->ui.lineEdit_gamma2->text();
	double D_new = text.toDouble();
	s->set_ptw(D_new); 
	//std::ostringstream os;os<<"Setting exponent for scalar field to: "<<D_new;report(os.str());
	//if(this->ui.comboBox_sf1->currentIndex()>0 & this->ui.comboBox_sf2->currentIndex()>0 & this->ui.comboBox_sf3->currentIndex()>0)
	//{
		this->s->update_scalar_fields();
		std::ostringstream os;
	os<<"\n"<<"Active region area fraction (Fold 1):"<<s->active_area_fraction_vis<<"\t("<<s->active_area_fraction<<")\n";
	os<<"\n"<<"Active region area fraction (Fold 2):"<<s->active_area_fraction_vis_2<<"\t("<<s->active_area_fraction_2<<")\n";
		report(os.str());		
		this->update_vtk();

	//}
	//this->update_vtk();
}
void SimGast::on_lineEdit_gamma3_editingFinished()
{
	const QString text  = this->ui.lineEdit_gamma3->text();
	double D_new = text.toDouble();
	s->set_phkb(D_new); 
	//std::ostringstream os;os<<"Setting exponent for scalar field to: "<<D_new;report(os.str());
		this->s->update_scalar_fields();
		std::ostringstream os;
	os<<"\n"<<"Active region area fraction (Fold 1):"<<s->active_area_fraction_vis<<"\t("<<s->active_area_fraction<<")\n";
	os<<"\n"<<"Active region area fraction (Fold 2):"<<s->active_area_fraction_vis_2<<"\t("<<s->active_area_fraction_2<<")\n";
		report(os.str());		
		this->update_vtk();
	//this->update_vtk();
}




void SimGast::on_lineEdit_stiff_fac_2_editingFinished()
{
	const QString text  = this->ui.lineEdit_stiff_fac_2->text();
	double D_new = text.toDouble();
	s->set_mp_fac_2(D_new); 
	//std::ostringstream os;os<<"Setting active region stiffness factor to "<<D_new;report(os.str());
	//this->update_vtk();
}
void SimGast::on_lineEdit_active_region_curvature_2_editingFinished()
{
	const QString text  = this->ui.lineEdit_active_region_curvature_2->text();
	double D_new = text.toDouble();
	s->set_mp_co_2(D_new); 
	//std::ostringstream os;os<<"Setting active region preferred curvature to "<<D_new;report(os.str());
	//this->update_vtk();
}


void SimGast::on_lineEdit_mp_cutoff_2_editingFinished()
{
	const QString text  = this->ui.lineEdit_mp_cutoff_2->text();
	double D_new = text.toDouble();
	this->s->set_mp_cutoff_2(D_new); 
	this->s->update_scalar_fields();
	std::ostringstream os;
	os<<"\n"<<"Active region area fraction (Fold 1):"<<s->active_area_fraction_vis<<"\t("<<s->active_area_fraction<<")\n";
	os<<"\n"<<"Active region area fraction (Fold 2):"<<s->active_area_fraction_vis_2<<"\t("<<s->active_area_fraction_2<<")\n";
	report(os.str());		
	this->update_vtk();
}
void SimGast::on_lineEdit_gamma1_2_editingFinished()
{
	const QString text  = this->ui.lineEdit_gamma1_2->text();
	double D_new = text.toDouble();
	s->set_psf1_2(D_new); 
		this->s->update_scalar_fields();
		std::ostringstream os;
		os<<"\n"<<"Active region area fraction (Fold 1):"<<s->active_area_fraction<<"\n";
		os<<"\n"<<"Active region area fraction (Fold 2):"<<s->active_area_fraction_2;
		report(os.str());		
		this->update_vtk();
}
void SimGast::on_lineEdit_gamma2_2_editingFinished()
{
	const QString text  = this->ui.lineEdit_gamma2_2->text();
	double D_new = text.toDouble();
	s->set_psf2_2(D_new); 
		this->s->update_scalar_fields();
		std::ostringstream os;
		os<<"\n"<<"Active region area fraction (Fold 1):"<<s->active_area_fraction<<"\n";
		os<<"\n"<<"Active region area fraction (Fold 2):"<<s->active_area_fraction_2;
		report(os.str());		
		this->update_vtk();
}
void SimGast::on_lineEdit_gamma3_2_editingFinished()
{
	const QString text  = this->ui.lineEdit_gamma3_2->text();
	double D_new = text.toDouble();
	s->set_psf3_2(D_new); 
		this->s->update_scalar_fields();
		std::ostringstream os;
		os<<"\n"<<"Active region area fraction (Fold 1):"<<s->active_area_fraction<<"\n";
		os<<"\n"<<"Active region area fraction (Fold 2):"<<s->active_area_fraction_2;
		report(os.str());		
		this->update_vtk();
}


void SimGast::update_comboBox_sf()
{
	if(s->srf_m->sf_tags.size())
	{
		this->ui.comboBox_sf1->blockSignals(true);
		this->ui.comboBox_sf1->clear();
		this->ui.comboBox_sf2->blockSignals(true);
		this->ui.comboBox_sf2->clear();
		this->ui.comboBox_sf3->blockSignals(true);
		this->ui.comboBox_sf3->clear();
		this->ui.comboBox_sf1_2->blockSignals(true);
		this->ui.comboBox_sf1_2->clear();
		this->ui.comboBox_sf2_2->blockSignals(true);
		this->ui.comboBox_sf2_2->clear();
		this->ui.comboBox_sf3_2->blockSignals(true);
		this->ui.comboBox_sf3_2->clear();
		
		this->ui.comboBox_sf1->addItem("none");
		this->ui.comboBox_sf2->addItem("none");
		this->ui.comboBox_sf3->addItem("none");
		this->ui.comboBox_sf1_2->addItem("none");
		this->ui.comboBox_sf2_2->addItem("none");
		this->ui.comboBox_sf3_2->addItem("none");
		
		for(int i=0;i<this->s->srf_m->sf_tags.size();i++)
		{
			this->ui.comboBox_sf1->addItem(s->srf_m->sf_tags[i].c_str());
			this->ui.comboBox_sf2->addItem(s->srf_m->sf_tags[i].c_str());
			this->ui.comboBox_sf3->addItem(s->srf_m->sf_tags[i].c_str());
			this->ui.comboBox_sf1_2->addItem(s->srf_m->sf_tags[i].c_str());
			this->ui.comboBox_sf2_2->addItem(s->srf_m->sf_tags[i].c_str());
			this->ui.comboBox_sf3_2->addItem(s->srf_m->sf_tags[i].c_str());
		}
		this->ui.comboBox_sf1->blockSignals(false);
		this->ui.comboBox_sf2->blockSignals(false);
		this->ui.comboBox_sf3->blockSignals(false);
		this->ui.comboBox_sf1_2->blockSignals(false);
		this->ui.comboBox_sf2_2->blockSignals(false);
		this->ui.comboBox_sf3_2->blockSignals(false);
	}
}
void SimGast::on_comboBox_sf1_currentIndexChanged(int val)
{
	this->s->sf1_vis.clear();
	s->sf1_vis.resize(gsmv->n_points);
	if(val>3)
	{
		gsmv->sfGen(s->srf_g->sfmx.col(val-4));	// generates the scalar field and stores it in gsmv->sf
		this->s->sfindx[0] = val-4;
	}
	else
	{
		this->s->sfindx[0] = -1;
		for(int i = 0;i<gsmv->n_points;i++){gsmv->sf[i] = 0;}
	}
	for(int i = 0;i<gsmv->n_points;i++){s->sf1_vis[i] = gsmv->sf[i];}	// copy contents of scalar field into scalar field visualization vector
	this->update_vtk();
}
void SimGast::on_comboBox_sf2_currentIndexChanged(int val)
{
	this->s->sf2_vis.clear();
	s->sf2_vis.resize(gsmv->n_points);
	if(val>3)
	{
		gsmv->sfGen(s->srf_g->sfmx.col(val-4));	// generates the surface and stores it in gsmv->sf
		this->s->sfindx[1] = val-4;
	}
	else
	{
		this->s->sfindx[1] = -1;
		for(int i = 0;i<gsmv->n_points;i++){gsmv->sf[i] = 0;}
	}
	// copy contents of scalar field into scalar field visualization vector
	for(int i = 0;i<gsmv->n_points;i++){s->sf2_vis[i] = gsmv->sf[i];}
	this->update_vtk();
}
void SimGast::on_comboBox_sf3_currentIndexChanged(int val)
{
	this->s->sf3_vis.clear();
	s->sf3_vis.resize(gsmv->n_points);
	if(val>3)
	{
		gsmv->sfGen(s->srf_g->sfmx.col(val-4));	// generates the surface and stores it in gsmv->sf
		this->s->sfindx[2] = val-4;
	}
	else
	{
		this->s->sfindx[2] = -1;
		for(int i = 0;i<gsmv->n_points;i++){gsmv->sf[i] = 0;}
	}
	// copy contents of scalar field into scalar field visualization vector
	s->sf3_vis.resize(gsmv->n_points);
	for(int i = 0;i<gsmv->n_points;i++){s->sf3_vis[i] = gsmv->sf[i];}
	this->update_vtk();
}


void SimGast::on_comboBox_sf1_2_currentIndexChanged(int val)
{
	this->s->sf1_2_vis.clear();
	s->sf1_2_vis.resize(gsmv->n_points);
	if(val>3)
	{
		gsmv->sfGen(s->srf_g->sfmx.col(val-4));	// generates the scalar field and stores it in gsmv->sf
		this->s->sfindx[3] = val-4;
	}
	else
	{
		this->s->sfindx[3] = -1;
		for(int i = 0;i<gsmv->n_points;i++){gsmv->sf[i] = 0;}
	}
	for(int i = 0;i<gsmv->n_points;i++){s->sf1_2_vis[i] = gsmv->sf[i];}	// copy contents of scalar field into scalar field visualization vector
	this->update_vtk();
}
void SimGast::on_comboBox_sf2_2_currentIndexChanged(int val)
{
	this->s->sf2_2_vis.clear();
	s->sf2_2_vis.resize(gsmv->n_points);
	if(val>3)
	{
		gsmv->sfGen(s->srf_g->sfmx.col(val-4));	// generates the surface and stores it in gsmv->sf
		this->s->sfindx[4] = val-4;
	}
	else
	{
		this->s->sfindx[4] = -1;
		for(int i = 0;i<gsmv->n_points;i++){gsmv->sf[i] = 0;}
	}
	// copy contents of scalar field into scalar field visualization vector
	s->sf2_2_vis.resize(gsmv->n_points);
	for(int i = 0;i<gsmv->n_points;i++){s->sf2_2_vis[i] = gsmv->sf[i];}
	this->update_vtk();
}
void SimGast::on_comboBox_sf3_2_currentIndexChanged(int val)
{
	this->s->sf3_2_vis.clear();
	s->sf3_2_vis.resize(gsmv->n_points);
	if(val>3)
	{
		gsmv->sfGen(s->srf_g->sfmx.col(val-4));	// generates the surface and stores it in gsmv->sf
		this->s->sfindx[5] = val-4;
	}
	else
	{
		this->s->sfindx[5] = -1;
		for(int i = 0;i<gsmv->n_points;i++){gsmv->sf[i] = 0;}
	}
	// copy contents of scalar field into scalar field visualization vector
	s->sf3_2_vis.resize(gsmv->n_points);
	for(int i = 0;i<gsmv->n_points;i++){s->sf3_2_vis[i] = gsmv->sf[i];}
	this->update_vtk();
}

void SimGast::update_all_sf()
{
	int val = 0;

	this->s->sf1_vis.clear();
	this->s->sf1_vis.resize(gsmv->n_points);
	val = this->s->sfindx[0] + 4;
	if(val>3)
	{
		gsmv->sfGen(s->srf_g->sfmx.col(val-4));	// generates the scalar field and stores it in gsmv->sf
		//this->s->sfindx[0] = val-4;
		for(int i = 0;i<gsmv->n_points;i++){s->sf1_vis[i] = gsmv->sf[i];}	// copy contents of scalar field into scalar field visualization vector
	}
	else{for(int i = 0;i<gsmv->n_points;i++){s->sf1_vis[i] = 0.0;}}
	
	this->s->sf2_vis.clear();
	s->sf2_vis.resize(gsmv->n_points);
	val = this->s->sfindx[1] + 4;
	if(val>3)
	{
		gsmv->sfGen(s->srf_g->sfmx.col(val-4));	// generates the surface and stores it in gsmv->sf
		//this->s->sfindx[1] = val-4;
		for(int i = 0;i<gsmv->n_points;i++){s->sf2_vis[i] = gsmv->sf[i];}// copy contents of scalar field into scalar field visualization vector
	}
	else{for(int i = 0;i<gsmv->n_points;i++){s->sf2_vis[i] = 0.0;}}

	this->s->sf3_vis.clear();
	s->sf3_vis.resize(gsmv->n_points);
	val = this->s->sfindx[2] + 4;
	if(val>3)
	{
		gsmv->sfGen(s->srf_g->sfmx.col(val-4));	// generates the surface and stores it in gsmv->sf
		//this->s->sfindx[2] = val-4;
		for(int i = 0;i<gsmv->n_points;i++){s->sf3_vis[i] = gsmv->sf[i];}	// copy contents of scalar field into scalar field visualization vector
	}
	else{for(int i = 0;i<gsmv->n_points;i++){s->sf3_vis[i] = 0.0;}}

	///////// now do the scalar fields relevant to the second fold
	this->s->sf1_2_vis.clear();
	this->s->sf1_2_vis.resize(gsmv->n_points);
	val = this->s->sfindx[3] + 4;
	if(val>3)
	{
		gsmv->sfGen(s->srf_g->sfmx.col(val-4));	// generates the scalar field and stores it in gsmv->sf
		//this->s->sfindx[0] = val-4;
	}
	for(int i = 0;i<gsmv->n_points;i++){s->sf1_2_vis[i] = gsmv->sf[i];}	// copy contents of scalar field into scalar field visualization vector
	
	this->s->sf2_2_vis.clear();
	this->s->sf2_2_vis.resize(gsmv->n_points);
	val = this->s->sfindx[4] + 4;
	if(val>3)
	{
		gsmv->sfGen(s->srf_g->sfmx.col(val-4));	// generates the surface and stores it in gsmv->sf
		//this->s->sfindx[1] = val-4;
	}
	// copy contents of scalar field into scalar field visualization vector
	for(int i = 0;i<gsmv->n_points;i++){s->sf2_2_vis[i] = gsmv->sf[i];}

	this->s->sf3_2_vis.clear();
	this->s->sf3_2_vis.resize(gsmv->n_points);
	val = this->s->sfindx[5] + 4;
	if(val>3)
	{
		gsmv->sfGen(s->srf_g->sfmx.col(val-4));	// generates the surface and stores it in gsmv->sf
		//this->s->sfindx[2] = val-4;
	}
	// copy contents of scalar field into scalar field visualization vector
	for(int i = 0;i<gsmv->n_points;i++){s->sf3_2_vis[i] = gsmv->sf[i];}

	////////////////////////////////////////////////////////////

	this->update_vtk();

}


// constraints and geometry
void SimGast::on_lineEdit_VMscale_editingFinished()
{
	const QString text  = this->ui.lineEdit_VMscale->text();
	double D_new = text.toDouble();
	s->set_VMscale(D_new); 
	//std::ostringstream os;os<<"Setting constraining surface scale \nrelative to outer shell to: "<<D_new;report(os.str());
	this->generate_vm_actor();
	this->update_vtk();
}
void SimGast::on_lineEdit_Lmax_editingFinished()
{
	const QString text  = this->ui.lineEdit_Lmax->text();
	int D_new = text.toInt();
	this->set_L_max(D_new);
}
void SimGast::on_lineEdit_gdim_editingFinished()
{
	const QString text  = this->ui.lineEdit_gdim->text();
	int D_new = text.toInt();
	if(D_new!=this->s->srf_m->b->dim)
	{
		std::ostringstream os;os<<"Setting Gaussian base point mesh to "<<D_new<<"  ..... ";report(os.str());
		this->s->set_new_gdim(D_new);
		this->gdim = D_new;
		os<<"Done!";report(os.str());
	}
}
void SimGast::on_lineEdit_gLmax_editingFinished()
{
	const QString text  = this->ui.lineEdit_gLmax->text();
	int D_new = text.toInt();
	if(D_new!=this->max_gL_max)
	{
		//s->set_gL_max(D_new); 
		std::ostringstream os;os<<"Setting maximum non-zero L to: "<<D_new;report(os.str());

		//if(D_new<s->get_gL_max())
		//if(D_new<this->max_gL_max & D_new>0)
		//{
			this->set_gLmax(D_new);
			//
			this->current_view_type = 1;		// set view type to display the scalar field on the midplane
			this->update_vtk();
		//}
		//else
		//if(D_new>this->max_gL_max | D_new<=0)
		//{
		//	ui.lineEdit_gLmax->setText(QString::number(double(this->max_gL_max)));
		//	string str = "L_max for scalar fields must be lower than current value and larger than zero.";
		//	message(str);
			
		//}
	}
}




//File menu, import and exports
void SimGast::on_actionOpen_session_triggered()
{
	QString filename = QFileDialog::getOpenFileName(0, "Open SimGast file ...",  "*.sgs");
	std::string filestr = (const char*)filename.toLatin1();
	if(!(filestr.length()==0)) 
	{
		ifstream in(filestr.c_str());
		if (!in)
		{	QMessageBox msgBox;
			msgBox.setText("Failed to open file for reading.");
			msgBox.setInformativeText("Press OK to continue.");
			msgBox.setStandardButtons(QMessageBox::Ok);
			msgBox.setDefaultButton(QMessageBox::Ok);
			int ret = msgBox.exec();
			return;
		}
		QProgressDialog progress("Importing and reconfiguring SimGast...", "Abort", 0, 100, 0);
		progress.setWindowModality(Qt::WindowModal);QApplication::processEvents();

		std::ostringstream os;os<<"\nImporting session from disc. Please wait ...."<<std::endl;report(os.str());os.str("");
		int version = 0;
		std::string tstr = "";

		in>>tstr;in >> version;
		os<<"\nFound version: "<<version<<"...."<<std::endl;report(os.str());os.str("");
		progress.setValue(10);
		progress.setLabelText("Importing shell object and generating basis. Please wait .....                   ");QApplication::processEvents();
		in>>tstr;in >> *s;
		// import all other application specific variables
		progress.setValue(50);
		progress.setLabelText("Reconfiguring application 1 of 2. Please wait .....                              ");QApplication::processEvents();
		in>>tstr;in>>this->opt_intersection_test_vm;
		in>>tstr;in>>this->opt_intersection_test_inner;
		in>>tstr;in>>this->opt_intersection_test_outer;
		in>>tstr;in>>this->show_hi_res;
		in>>tstr;in>>this->opt_plot_step;

		if( version>1000)
		{
			bool check;
			this->ui.checkBox_C_o_model->blockSignals(true);
			this->ui.checkBox_symmetry->blockSignals(true);
			in>>tstr;in>>check;this->ui.checkBox_C_o_model->setChecked(check);
			in>>tstr;in>>check;this->ui.checkBox_symmetry->setChecked(check);
			this->ui.checkBox_C_o_model->blockSignals(false);
			this->ui.checkBox_symmetry->blockSignals(false);
		}
		int tw_indx;
		in>>tw_indx;	// for the optimization tab widget
		
		int val;
		double dval;
		////// SUBPLEX optimization
		in>>val;ui.lineEdit_SBPLX_Lmax_fit->setText(QString::number(val));
		in>>val;ui.lineEdit_SBPLX_max_fun->setText(QString::number(val));
		in>>dval;ui.lineEdit_SBPLX_gamma_inter->setText(QString::number(dval));
		in>>dval;ui.lineEdit_SBPLX_dx->setText(QString::number(dval));
		in>>dval;ui.lineEdit_SBPLX_gamma_v->setText(QString::number(dval));
		////// COBYLA optimization
		in>>val;ui.lineEdit_COBYLA_Lmax_fit->setText(QString::number(val));
		in>>val;ui.lineEdit_COBYLA_max_fun->setText(QString::number(val));
		in>>dval;ui.lineEdit_COBYLA_xtol->setText(QString::number(dval));
		in>>dval;ui.lineEdit_COBYLA_dx->setText(QString::number(dval));
		in>>dval;ui.lineEdit_COBYLA_vtol->setText(QString::number(dval));
		/////// ISRES optimization
		in>>val;ui.lineEdit_ISRES_Lmax_fit->setText(QString::number(val));
		in>>val;ui.lineEdit_ISRES_max_fun->setText(QString::number(val));
		in>>val;ui.lineEdit_ISRES_pop->setText(QString::number(val));
		in>>dval;ui.lineEdit_ISRES_lb->setText(QString::number(dval));
		in>>dval;ui.lineEdit_ISRES_ub->setText(QString::number(dval));
		in>>dval;ui.lineEdit_ISRES_vtol->setText(QString::number(dval));

		in.close();
		progress.setValue(80);
		progress.setLabelText("Importing shell object 2 of 2. Please wait .....                                 ");QApplication::processEvents();
		// update state of the shell object
		s->set_L_max(this->Lmax);
		//s->srf_g->update();
		this->s->srf_m->center_to_zero();
		this->s->srf_u->center_to_zero();
		this->s->srf_g->center_to_zero();

		//this->s->set_undeformed_shape();		// use srf_m as the undeformed shape srf_u
		this->update_comboBox_sf();
		s->srf_m->needsUpdating = YES;
		s->srf_u->needsUpdating = YES;
		s->fix_normals_sign_known =0;
		s->fix_normals_sign = 1;
		this->s->srf_m->update();
		this->s->srf_m->update_tri();
		this->s->srf_u->update();
		this->s->srf_u->update_tri();
		this->s->srf_g->update();
		this->s->srf_g->update_tri();


		/// initialize the visualization mesh objects
		/**/
		progress.setValue(90);
		progress.setLabelText("Updating visualization components. Please wait .....                                 ");QApplication::processEvents();
		spherical_mesh* tmp_smv = this->smv;
		spherical_mesh* tmp_gsmv = this->gsmv;
		smv = new spherical_mesh(this->s->srf_m->b->L_max, tmp_smv->tri_n);
		gsmv= new spherical_mesh(this->s->srf_g->b->L_max, tmp_gsmv->tri_n);
		delete tmp_smv;
		delete tmp_gsmv;
		
		s->sf1_vis.resize(gsmv->n_points);
		s->sf2_vis.resize(gsmv->n_points);
		s->sf3_vis.resize(gsmv->n_points);
		generate_vm_actor();
		cscale = 1.0/cscale_calc(this->s->srf_m);
		_axes_actor1->SetTotalLength(_scale_fac/cscale, _scale_fac/cscale, _scale_fac/cscale);
		// update the state of the user interface according to the values read from this file
		update_comboBox_sf();		// fills the pulldown menus of the scalar fields
		this->ui.comboBox_sf1->blockSignals(true);
		this->ui.comboBox_sf2->blockSignals(true);
		this->ui.comboBox_sf3->blockSignals(true);
		this->ui.comboBox_sf1_2->blockSignals(true);
		this->ui.comboBox_sf2_2->blockSignals(true);
		this->ui.comboBox_sf3_2->blockSignals(true);
		ui.comboBox_sf1->setCurrentIndex(s->sfindx[0] + 4);
		ui.comboBox_sf2->setCurrentIndex(s->sfindx[1] + 4);
		ui.comboBox_sf3->setCurrentIndex(s->sfindx[2] + 4);
		ui.comboBox_sf1_2->setCurrentIndex(s->sfindx[3] + 4);
		ui.comboBox_sf2_2->setCurrentIndex(s->sfindx[4] + 4);
		ui.comboBox_sf3_2->setCurrentIndex(s->sfindx[5] + 4);
		this->ui.comboBox_sf1->blockSignals(false);
		this->ui.comboBox_sf2->blockSignals(false);
		this->ui.comboBox_sf3->blockSignals(false);
		this->ui.comboBox_sf1_2->blockSignals(false);
		this->ui.comboBox_sf2_2->blockSignals(false);
		this->ui.comboBox_sf3_2->blockSignals(false);

		ui.checkBox_vm->setChecked(opt_intersection_test_vm);
		ui.checkBox_inner->setChecked(opt_intersection_test_inner);
		ui.checkBox_outer->setChecked(opt_intersection_test_outer);
		ui.checkBox_HiRes->setChecked(show_hi_res);
		ui.lineEdit_opt_plot_frequency->setText(QString::number(opt_plot_step));
		ui.tabWidget->setCurrentIndex(tw_indx);
		//ui.checkBox_symmetry->setChecked(0);
		this->on_checkBox_C_o_model_stateChanged(this->ui.checkBox_C_o_model->checkState()); // to actually do the function call
		this->on_checkBox_symmetry_stateChanged(this->ui.checkBox_symmetry->checkState()); // to trigger the function call

		
		this->s->update_scalar_fields();
		update_all_sf();
		update_vtk();
		os<<"\nDone!!"<<std::endl;report(os.str());os.str("");
		progress.setValue(100);
	}
	/*
	QMessageBox msgBox;
	msgBox.setText("The function has not been implemented yet.");
	msgBox.setInformativeText("Press OK to continue.");
	msgBox.setStandardButtons(QMessageBox::Ok);
	msgBox.setDefaultButton(QMessageBox::Ok);
	int ret = msgBox.exec();
	*/
}
void SimGast::on_actionSave_session_triggered()
{
	QString filename = QFileDialog::getSaveFileName(0, "Save session to file ...",  "*.sgs");
	std::string filestr = (const char*)filename.toLatin1();
	if(!(filestr.length()==0)) 
	{
		ofstream out(filestr.c_str());
		if (!out)
		{	QMessageBox msgBox;
			msgBox.setText("Failed to open file for writing.");
			msgBox.setInformativeText("Press OK to continue.");
			msgBox.setStandardButtons(QMessageBox::Ok);
			msgBox.setDefaultButton(QMessageBox::Ok);
			int ret = msgBox.exec();
			return;
		}
		out<<"version"<<std::endl;out<<this->version<<std::endl;
		out<<"shell_object>>>>>>>>>>>>>>>>>>>>>>>>>>>>"<<std::endl;
			out<<*s;		// saves the shell object and all its members

			/// now export all other application specific variables

			out<<"intersection_test_vm"<<std::endl;out<<this->ui.checkBox_vm->isChecked()<<std::endl;//out<<this->opt_intersection_test_vm <<std::endl;
			out<<"intersection_test_inner"<<std::endl;out<<this->ui.checkBox_inner->isChecked()<<std::endl;//out<<this->opt_intersection_test_inner<<std::endl;
			out<<"intersection_test_outer"<<std::endl;out<<this->ui.checkBox_outer->isChecked()<<std::endl;//out<<this->opt_intersection_test_outer<<std::endl;
			out<<"show_hi_resolution"<<std::endl;out<<this->ui.checkBox_HiRes->isChecked()<<std::endl;//out<<this->show_hi_res<<std::endl;
			out<<"fitting_plot_step"<<std::endl;out<<this->opt_plot_step<<std::endl;
			out<<"model_for_C_o"<<std::endl;out<<this->ui.checkBox_C_o_model->isChecked()<<std::endl; //
			out<<"use_mirror_symmetry"<<std::endl;out<<this->ui.checkBox_symmetry->isChecked()<<std::endl; //

			out<<ui.tabWidget->currentIndex()<<std::endl;
			
			QString text;
			////// SUBPLEX parameters
			
			text  = this->ui.lineEdit_SBPLX_Lmax_fit->text();out<<text.toInt()<<std::endl;
			text  = this->ui.lineEdit_SBPLX_max_fun->text();out<<text.toInt()<<std::endl;
			text  = this->ui.lineEdit_SBPLX_gamma_inter->text();out<<text.toDouble()<<std::endl;
			text  = this->ui.lineEdit_SBPLX_dx->text();out<<text.toDouble()<<std::endl;
			text  = this->ui.lineEdit_SBPLX_gamma_v->text();out<<text.toDouble()<<std::endl;

			////// COBYLA parameters
			text  = this->ui.lineEdit_COBYLA_Lmax_fit->text();out<<text.toInt()<<std::endl;
			text  = this->ui.lineEdit_COBYLA_max_fun->text();out<<text.toInt()<<std::endl;
			text  = this->ui.lineEdit_COBYLA_xtol->text();out<<text.toDouble()<<std::endl;
			text  = this->ui.lineEdit_COBYLA_dx->text();out<<text.toDouble()<<std::endl;
			text  = this->ui.lineEdit_COBYLA_vtol->text();out<<text.toDouble()<<std::endl;
			
			////// ISRES parameters
			text  = this->ui.lineEdit_ISRES_Lmax_fit->text();out<<text.toInt()<<std::endl;
			text  = this->ui.lineEdit_ISRES_max_fun->text();out<<text.toInt()<<std::endl;
			text  = this->ui.lineEdit_ISRES_pop->text();out<<text.toInt()<<std::endl;
			text  = this->ui.lineEdit_ISRES_lb->text();out<<text.toDouble()<<std::endl;
			text  = this->ui.lineEdit_ISRES_ub->text();out<<text.toDouble()<<std::endl;
			text  = this->ui.lineEdit_ISRES_vtol->text();out<<text.toDouble()<<std::endl;

		/// close the file
		out.close();
		
	}
	/*
	QMessageBox msgBox;
	msgBox.setText("The function has not been implemented yet.");
	msgBox.setInformativeText("Press OK to continue.");
	msgBox.setStandardButtons(QMessageBox::Ok);
	msgBox.setDefaultButton(QMessageBox::Ok);
	int ret = msgBox.exec();
	*/
}
void SimGast::on_actionExit_application_triggered() {
  qApp->exit();
}
void SimGast::on_actionSave_shp3_triggered()
{
	QString filename = QFileDialog::getSaveFileName(0, "Save shape to file ...",  "*.shp3");
	std::string filestr = (const char*)filename.toLatin1();
	if(!(filestr.length()==0)) 
	{
		if(s->srf_m->b->L_max<s->srf_g->b->L_max)
		{
			this->set_L_max(s->srf_g->b->L_max);
		}
		else
		if(s->srf_m->b->L_max>s->srf_g->b->L_max)
		{
			this->set_gLmax(s->srf_m->b->L_max);
		}
		MatrixXd *tmp = &s->srf_m->sfmx;
		MatrixXd tmp2 = s->srf_g->sfmx;
		//tmp2.resize(s->srf_m->sfmx.rows(), s->srf_m->sfmx.cols());
		s->srf_m->sfmx = tmp2;
		s->srf_m->write(filestr.c_str());
		tmp->resize(0,0);
	}
}
void SimGast::on_actionExport_obj_triggered()
{
    /*
	QString filename = QFileDialog::getSaveFileName(0, "",  "*.obj");
	std::string filestr = (const char*)filename.toLatin1();
	if(!(filestr.length()==0))
	{
		int tmp = 0;
		if(current_view_type==0)	// i.e. if we are looking at a cut-shell
		{
			// generate the cut-shell
			vtkPolyData*		polydata_in	= vtkPolyData::New();
			vtkPolyData*		polydata_out= vtkPolyData::New();
			vtkPolyData*		polydata_fill= vtkPolyData::New();
			shell_cut2polydata1(this->s, polydata_in, polydata_out, polydata_fill);
			Mapper1_in->SetInputData(polydata_in);
			Mapper1_in->SetScalarRange(_min_sf_1, _max_sf_1);
			Actor1_in->SetMapper( Mapper1_in);
			Mapper1_out->SetInputData(polydata_out);
			Mapper1_out->SetScalarRange(_min_sf_1, _max_sf_1);
			Actor1_out->SetMapper( Mapper1_out);
			Mapper1_fill->SetInputData(polydata_fill);
			Mapper1_fill->SetScalarRange(_min_sf_1, _max_sf_1);
			Actor1_fill->SetMapper( Mapper1_fill);
			vtkExtractEdges *edges = vtkExtractEdges::New();
			edges->SetInputData(polydata_fill);
			vtkPolyDataMapper *edge_mapper = vtkPolyDataMapper::New();
			edge_mapper->SetInputData(edges->GetOutput());
			vtkActor *edge_actor = vtkActor::New();
			edge_actor->SetMapper(edge_mapper);
			edge_actor->GetProperty()->SetColor(1,0.5,0.5);
			// prepare the scene by removing the axis actor if present
			_ren1->RemoveAllViewProps(); // remove the existing actors
			if(axis_on_1)
			{
				_ren1->AddActor(_axes_actor1);		//  add the axis actor
			}
							
			_ren1->AddActor(this->Actor1_in);
			_ren1->AddActor(this->Actor1_out);
			_ren1->AddActor(this->Actor1_fill);//_ren1->AddActor(edge_actor);
			_ren1->ResetCamera();
			_ren1->GetRenderWindow()->Render();

		}
		else
			//(current_view_type==1)	// i.e. if we are looking at a cut-shell
		{
			//vtkPolyData*		polydata	= vtkPolyData::New();
			this->s->update();
			shell2polydata1();
			Mapper1_in->SetInputData(vis_polydata);
			//Mapper1_in->SetScalarRange(_min_sf_1, _max_sf_1);
			Actor1_in->SetMapper( Mapper1_in);
			_ren1->RemoveAllViewProps();
			if(this->axis_on_1){_ren1->AddActor(this->_axes_actor1);}
			_ren1->AddActor(this->Actor1_in);
			_ren1->GetRenderWindow()->Render();
		
		}
		// code below will be executed in any case
		vtkOBJExporter *obj = vtkOBJExporter::New();
		obj->SetInput(_ren1->GetRenderWindow());
		obj->SetFilePrefix((filename.toStdString()).c_str());
		obj->Write();
	}
     */
}

void SimGast::on_pushButton_load_undeformed_clicked()
{
	// load the shape
	QString filename = QFileDialog::getOpenFileName(ui.vtkWidget1, "Open undeformed shape file ...",  "*.shp3");
	this->filestr = (const char*) filename.toLatin1();//filename.toStdString();
	if(!(filestr.length()==0)){	this->import_shape_from_disc(filestr);}
	ui.checkBox_symmetry->setChecked(0);

}
void SimGast::on_pushButton_load_starting_clicked()
{
	// load the shape
	QString filename = QFileDialog::getOpenFileName(ui.vtkWidget1, "Open starting shape file ...",  "*.shp3");
	this->filestr = (const char*) filename.toLatin1();//filename.toStdString();
	if(!(filestr.length()==0))
	{
		this->s->srf_m->read_trunc(this->filestr.c_str());		// read_trunc causes truncation of series or filling with zeros to keep L_max the same
		this->s->srf_m->center_to_zero();
		s->srf_m->needsUpdating = YES;
		s->srf_m->update();
		s->fix_normals_sign_known =0;
		s->fix_normals_sign = 1;
		s->update();
		cscale = 1.0/cscale_calc(this->s->srf_m);
		_axes_actor1->SetTotalLength(_scale_fac/cscale, _scale_fac/cscale, _scale_fac/cscale);
		update_vtk();
	}
}
QString SimGast::mimeData_2_fileName( const QMimeData * mimeData )
{
	if ( mimeData->hasUrls() )
	{
		foreach ( const QUrl & url, mimeData->urls() )
		{
			QString str = url.toLocalFile();
			if ( str.isEmpty() == false )
			{
				return str;
			}
		}
	}
	return QString();
}
void SimGast::dragEnterEvent( QDragEnterEvent * event )
{
	if ( mimeData_2_fileName( event->mimeData() ).isEmpty() == false )
	{
		event->acceptProposedAction();
	}
}
void SimGast::dropEvent( QDropEvent * event )
{
	QString filename = mimeData_2_fileName( event->mimeData() );
	this->filestr = (const char*) filename.toLatin1();//filename.toStdString();
	if ( filename.isEmpty() == false ){this->import_shape_from_disc(filestr);}
}

// reporting and visualization
void SimGast::on_pushButton_status_clicked()
{
	std::ostringstream os;s->disp(os,*s);//report(os.str());
	if(this->ui.comboBox_sf1->currentIndex()>0 & this->ui.comboBox_sf2->currentIndex()>0 & this->ui.comboBox_sf3->currentIndex()>0)
	{
		int old_tri_n = this->s->srf_m->sm->tri_n;
		this->s->srf_m->set_new_spherical_mesh(tri_n_intersection_test);		// we use a low resolution mesh for intersection tests
		double E = this->s->get_energy(1);
		this->s->srf_m->set_new_spherical_mesh(old_tri_n);
		//this->s->srf_m->update();
		os<<"\n"<<"Active region area fraction:"<<s->active_area_fraction<<std::endl;
		//os<<"Vout Vin Aout Ain:"<<<<"\n"<<s->vin<<"\t"<<s->aout<<"\t"<<s->ain<<"\t"<<std::endl;
		os<<"Volume enclosed by outer surface\t"<<s->vout<<std::endl;
		os<<"Area of outer surface\t"<<s->aout<<std::endl;
		os<<"Volume enclosed by inner surface\t\t"<<s->vin<<std::endl;
		os<<"Area of inner surface\t"<<s->ain<<std::endl;
		os<<"Shear and stretch: "<<this->s->E_nHook<<"\tBending: "<<this->s->Eb<<"\nTotal energy:\n"<<E<<std::endl;

	}
	report(os.str());

}
void SimGast::on_pushButton_clear_clicked(){this->ui.textBrowser_report->clear();}
void SimGast::on_pushButton_rot_x_clicked()
{
	double ang[3];
	ang[0] = PI/2;
	ang[1] = PI/2;
	ang[2] = -PI/2;
	s->srf_m->rotate_around_self(ang);	s->srf_m->needsUpdating = YES;
	//s->srf_u->rotate_around_self(ang);	s->srf_u->needsUpdating = YES;this->generate_vm_actor();
	std::ostringstream os;os<<"Rotating current shape only around x-axis by pi/2 ";report(os.str());
	
	this->update_vtk();
}
void SimGast::on_pushButton_rot_y_clicked()
{
	double ang[3];
	ang[0] = PI/2;
	ang[1] = 0.0;
	ang[2] = 0.0;
	s->srf_m->rotate_around_self(ang);
	s->srf_m->needsUpdating = YES;
	//s->srf_u->rotate_around_self(ang);	s->srf_u->needsUpdating = YES;this->generate_vm_actor();
	std::ostringstream os;os<<"Rotating current shape only around y-axis by pi/2 ";report(os.str());
	this->update_vtk();
}
void SimGast::on_pushButton_rot_z_clicked()
{
	double ang[3];
	ang[0] = 0.0;
	ang[1] = PI/2;
	ang[2] = 0.0;
	s->srf_m->rotate_around_self(ang);	s->srf_m->needsUpdating = YES;
	//s->srf_u->rotate_around_self(ang);	s->srf_u->needsUpdating = YES;this->generate_vm_actor();
	std::ostringstream os;os<<"Rotating current shape only around z-axis by pi/2 ";report(os.str());
	this->update_vtk();

	
}
void SimGast::on_pushButton_set_undef_clicked()
{

	std::ostringstream os;os<<"Setting the undeformed shape to current one";report(os.str());
	this->s->set_undeformed_shape();
	s->srf_u->needsUpdating = YES;
	this->s->srf_u->update();
	this->s->srf_u->update_tri();
	this->generate_vm_actor();
	this->update_vtk();

}

void SimGast::on_pushButton_reset_view_clicked()
{
	this->update_vtk();
	_ren1->ResetCamera();
	_ren1->GetRenderWindow()->Render();
}
void SimGast::on_pushButton_section_clicked()
{
	current_view_type = 0;
	this->update_vtk();
}
void SimGast::on_pushButton_curvature_clicked()
{
	current_view_type = 3;
	this->update_vtk();
}
void SimGast::on_pushButton_mid_plane_clicked()
{
	//double v1, v2, v3, p1, p2, p3;
	//this->camera1->GetViewUp(v1,v2,v3);
	//this->camera1->GetPosition(p1,p2,p3);
	//std::ostringstream os;os<<"\n";
	//os<<"Camera angle"<<this->camera1->GetViewAngle()<<"\n";
	//os<<"Camera view up"<<v1<<"\t"<<v2<<"\t"<<v3<<"\n";
	//os<<"Camera position"<<p1<<"\t"<<p2<<"\t"<<p3<<"\n";
	//report(os.str());
	current_view_type = 1;
	this->update_vtk();
}
void SimGast::on_pushButton_axis_toggle_clicked()
{
	if(this->axis_on_1){this->axis_on_1 = 0;}
	else
	{this->axis_on_1 = 1;}
	this->update_vtk();
}
void SimGast::on_pushButton_cube_axis_toggle_clicked()
{
	if(this->cube_axis_actor1_on){this->cube_axis_actor1_on = 0;}
	else
	{this->cube_axis_actor1_on = 1;}
	this->update_vtk();
}

void SimGast::on_pushButton_active_region_view_clicked()
{
	//if(this->ui.comboBox_sf1->currentIndex()>0 & this->ui.comboBox_sf2->currentIndex()>0 & this->ui.comboBox_sf3->currentIndex()>0)
	//{
	if (this->s->sfindx[0] >-1 | this->s->sfindx[1] >-1 | this->s->sfindx[2] >-1 | this->s->sfindx[3] >-1 | this->s->sfindx[4] >-1 | this->s->sfindx[5] >-1)
	{
		if(current_view_type !=2)
		{
			current_view_type = 2;
		}
		else
		{
			current_view_type = 1;
		}
		this->s->update_scalar_fields();
		std::ostringstream os;
		os<<"\n"<<"Active region area fraction (Fold 1):"<<s->active_area_fraction<<"\n";
		os<<"\n"<<"Active region area fraction (Fold 2):"<<s->active_area_fraction_2;
		report(os.str());		
		this->update_vtk();
	}
	//}
	//else
	//{
	//	this->message("Please select three scalar fields first.");
	//}
}

void SimGast::on_pushButton_constraint_surface_clicked()
{
	
	if(vm_actor_on)
	{
		vm_actor_on = false;
		vmActor->SetVisibility(false);
	}
	else
	{
		vm_actor_on = true;
		vmActor->SetVisibility(true);
	}
	this->update_vtk();
}
// general settings
void SimGast::on_lineEdit_opt_plot_frequency_editingFinished()
{
	const QString text  = this->ui.lineEdit_opt_plot_frequency->text();
	this->opt_plot_step = text.toInt();
}





void SimGast::on_checkBox_vm_stateChanged(int state)
{
	if(state==2){opt_intersection_test_vm = true;}
	else{opt_intersection_test_vm = false;}
}
void SimGast::on_checkBox_inner_stateChanged(int state)
{
	if(state==2){opt_intersection_test_inner = true;}
	else{opt_intersection_test_inner = false;}
}
void SimGast::on_checkBox_outer_stateChanged(int state)
{
	if(state==2){opt_intersection_test_outer = true;}
	else{opt_intersection_test_outer = false;}
}
void SimGast::on_checkBox_HiRes_stateChanged(int state)
{
	if(state==2){show_hi_res = true;}
	else{show_hi_res = false;}
}
void SimGast::on_checkBox_symmetry_stateChanged(int state)
{	std::ostringstream os;
	if(state==2)
	{
		this->enforce_symmetry();
		os<<"Mirror symmetry plane enforcement is in effect.\nSome coefficients have been affected.\nOpimization routine considers reduced set only.\n"<<std::endl;report(os.str());os.str("");
		update_vtk();
	}
	else
	{
		os<<"Note: opimization will consider the full coefficient set.\n"<<std::endl;report(os.str());os.str("");
		int nc = this->s->srf_m->xc.rows();
		ixc.resize(nc);std::fill(ixc.begin(), ixc.end(), true);
		iyc.resize(nc);std::fill(iyc.begin(), iyc.end(), true);
		izc.resize(nc);std::fill(izc.begin(), izc.end(), true);
		this->ui.checkBox_symmetry_vfi->blockSignals(true);
		this->ui.checkBox_symmetry_vfi->setChecked(false);
		this->ui.checkBox_symmetry_vfi->blockSignals(false);

	}
}
void SimGast::on_checkBox_symmetry_vfi_stateChanged(int state)
{	std::ostringstream os;
/*	if(state==2)
	{
		this->enforce_symmetry_vfi();
		os<<"VFI forming coefficients basis pruning.\nMirror symmetry plane enforcement is also in effect.\nSome coefficients have been affected.\nOpimization routine considers reduced set only.\n"<<std::endl;report(os.str());os.str("");
		this->ui.checkBox_symmetry->blockSignals(true);
		this->ui.checkBox_symmetry->setChecked(true);
		this->ui.checkBox_symmetry->blockSignals(false);
		update_vtk();
	}
	else
	{
		os<<"Note: opimization will consider the full coefficient set.\n"<<std::endl;report(os.str());os.str("");
		int nc = this->s->srf_m->xc.rows();
		ixc.resize(nc);std::fill(ixc.begin(), ixc.end(), true);
		iyc.resize(nc);std::fill(iyc.begin(), iyc.end(), true);
		izc.resize(nc);std::fill(izc.begin(), izc.end(), true);
		this->ui.checkBox_symmetry->blockSignals(true);
		this->ui.checkBox_symmetry->setChecked(false);
		this->ui.checkBox_symmetry->blockSignals(false);
		
	}

*/
}

void SimGast::on_checkBox_C_o_model_stateChanged(int state)
{	std::ostringstream os;
	if(state==2)
	{
		this->s->set_C_o_model(1);
		os<<"Preferred curvature is set to that of the undeformed geometry\n"<<std::endl;report(os.str());os.str("");
	}
	else
	{
		this->s->set_C_o_model(0);
		os<<"Preferred curvature is set to zero (flat)\n"<<std::endl;report(os.str());os.str("");

	}
}

//////////////////////// Surface generation ////////////////
void SimGast::shell2polydata1()
{
	/*vtkCellArray *_polys2 = vtkCellArray::New();
	for (int i=0; i<this->s->srf_m->sm->n_faces; i++)
	{
		_polys2->InsertNextCell(3);
		_polys2->InsertCellPoint(this->s->srf_m->sm->f0[i]-1);
		_polys2->InsertCellPoint(this->s->srf_m->sm->f1[i]-1);
		_polys2->InsertCellPoint(this->s->srf_m->sm->f2[i]-1);
	}
	vtkPoints* points = vtkPoints::New();
	points->SetDataTypeToDouble();
	points->SetNumberOfPoints(s->srf_m->sm->n_points);
	for (int i=0; i<s->srf_g->sm->n_points; i++){	points->SetPoint(i, s->srf_m->sm->X[i][0], s->srf_m->sm->X[i][1], s->srf_m->sm->X[i][2]);}
	polydata->SetPoints(points);points->Delete();// Assign points and cells
	polydata->SetPolys(_polys2);//polys->Delete();
	*/
	this->smv->surfaceGenX(this->s->srf_m->xc, this->s->srf_m->yc, this->s->srf_m->zc); // use fast surface generation for X only
	vtkPoints* points = vtkPoints::New();
	points->SetDataTypeToDouble();
	points->SetNumberOfPoints(this->smv->n_points);
	for (int i=0; i<this->smv->n_points; i++)
	{	
		points->SetPoint(i, this->smv->Xu(i,0), this->smv->Xu(i,1), this->smv->Xu(i,2));/// fill in the points
	}
	vis_polydata->SetPoints(points);points->Delete();// Assign points
	////////////////////////////////////////////////
	// add the scalar fields if needed
	int nsize = s->sf3_vis.size();
	if(this->ui.comboBox_sf1->currentIndex()>3 || this->ui.comboBox_sf2->currentIndex()>3 || this->ui.comboBox_sf3->currentIndex()>3 ||this->ui.comboBox_sf1_2->currentIndex()>3 || this->ui.comboBox_sf2_2->currentIndex()>3 || this->ui.comboBox_sf3_2->currentIndex()>3)
	{
		_min_sf_1 = 0.0;
		_max_sf_1 = 0.0;
		//int nsize = 0;
		if(this->ui.comboBox_sf1->currentIndex()>3)
		{
			float tmp;
			for (int i=0; i<s->sf1_vis.size(); i++)
			{
 				if((this->s->sf1_vis[i])>_max_sf_1){_max_sf_1 = this->s->sf1_vis[i];}
				if((this->s->sf1_vis[i])<_min_sf_1){_min_sf_1 = this->s->sf1_vis[i];}
			}
			nsize = s->sf1_vis.size();
		}
		if(this->ui.comboBox_sf2->currentIndex()>3)
		{
			float tmp;
			for (int i=0; i<s->sf2_vis.size(); i++)
			{
 				if((this->s->sf2_vis[i])>_max_sf_1){_max_sf_1 = this->s->sf2_vis[i];}
				if((this->s->sf2_vis[i])<_min_sf_1){_min_sf_1 = this->s->sf2_vis[i];}
			}
			nsize = s->sf2_vis.size();
		}
		if(this->ui.comboBox_sf3->currentIndex()>3)
		{
			float tmp;
			for (int i=0; i<s->sf3_vis.size(); i++)
			{
 				if((this->s->sf3_vis[i])>_max_sf_1){_max_sf_1 = this->s->sf3_vis[i];}
				if((this->s->sf3_vis[i])<_min_sf_1){_min_sf_1 = this->s->sf3_vis[i];}
			}
			nsize = s->sf3_vis.size();
		}

		if(this->ui.comboBox_sf1_2->currentIndex()>3)
		{
			float tmp;
			for (int i=0; i<s->sf1_2_vis.size(); i++)
			{
 				if((this->s->sf1_2_vis[i])>_max_sf_1){_max_sf_1 = this->s->sf1_2_vis[i];}
				if((this->s->sf1_2_vis[i])<_min_sf_1){_min_sf_1 = this->s->sf1_2_vis[i];}
			}
			nsize = s->sf1_2_vis.size();
		}
		if(this->ui.comboBox_sf2_2->currentIndex()>3)
		{
			float tmp;
			for (int i=0; i<s->sf2_2_vis.size(); i++)
			{
 				if((this->s->sf2_2_vis[i])>_max_sf_1){_max_sf_1 = this->s->sf2_2_vis[i];}
				if((this->s->sf2_2_vis[i])<_min_sf_1){_min_sf_1 = this->s->sf2_2_vis[i];}
			}
			nsize = s->sf2_2_vis.size();
		}
		if(this->ui.comboBox_sf3_2->currentIndex()>3)
		{
			float tmp;
			for (int i=0; i<s->sf3_2_vis.size(); i++)
			{
 				if((this->s->sf3_2_vis[i])>_max_sf_1){_max_sf_1 = this->s->sf3_2_vis[i];}
				if((this->s->sf3_2_vis[i])<_min_sf_1){_min_sf_1 = this->s->sf3_2_vis[i];}
			}
			nsize = s->sf3_2_vis.size();
		}

		vtkUnsignedCharArray *scalars = vtkUnsignedCharArray::New();
		scalars->SetNumberOfComponents(3);
		scalars->SetNumberOfTuples(nsize  );	//assuming all scalar fields to hold the same number of elements
		//double fac = 255.0/_max_sf_1;
		double fac = 255.0/(_max_sf_1-_min_sf_1);
		//std::ostringstream os;
		//os<<"_max_sf_1: "<<_max_sf_1<<" --- _min_sf1: "<<_min_sf_1<<std::endl;
		//os<<"No. of points: "<<nsize<<"\n";
		for (int i=0; i<nsize; i++)
		{
			double s1 = 0.0;double s2 = 0.0; double s3 = 0.0;double s1_2 = 0.0;double s2_2 = 0.0; double s3_2 = 0.0;
			if(this->ui.comboBox_sf1->currentIndex()>3){s1 = (this->s->sf1_vis[i] - _min_sf_1)*fac;}
			if(this->ui.comboBox_sf2->currentIndex()>3){s2 = (this->s->sf2_vis[i] - _min_sf_1)*fac;}
			if(this->ui.comboBox_sf3->currentIndex()>3){s3 = (this->s->sf3_vis[i] - _min_sf_1)*fac;}
			
			if(this->ui.comboBox_sf1_2->currentIndex()>3){s1_2 = (this->s->sf1_2_vis[i] - _min_sf_1)*fac;}
			if(this->ui.comboBox_sf2_2->currentIndex()>3){s2_2 = (this->s->sf2_2_vis[i] - _min_sf_1)*fac;}
			if(this->ui.comboBox_sf3_2->currentIndex()>3){s3_2 = (this->s->sf3_2_vis[i] - _min_sf_1)*fac;}

			scalars->InsertTuple3(i, (s1 + s1_2)/2.0, (s2 + s2_2)/2.0, (s3 + s3_2)/2.0);
			//if(s1>=254.0){os<<i<<"\t"<<s1<<"\t"<<this->s->sf1_vis[i]<<"\t"<<s3<<"\n";}
		}
		//report(os.str());
		vis_polydata->GetPointData()->SetScalars(scalars);//scalars->Delete();//uncomment if there are curvature values associated with polydata
	}
	else
	{
		vtkUnsignedCharArray *scalars = vtkUnsignedCharArray::New();
		scalars->SetNumberOfComponents(3);
		scalars->SetNumberOfTuples(nsize  );	//assuming all scalar fields to hold the same number of elements
		for (int i=0; i<nsize; i++){scalars->InsertTuple3(i,255, 255, 255);}
		vis_polydata->GetPointData()->SetScalars(scalars);
	}

}
void SimGast::shell2polydata1_active_region()
{
		//std::ostringstream os;
		//os<<"min and max field values found (before normalization): "<<_min_sf_1<<"\t"<<_max_sf_1<<"\n";
		//os<<"No. of points MP_vis: "<<s->MP_vis.size()<<"\n";report(os.str());os.str("");
		vtkUnsignedCharArray *scalars = vtkUnsignedCharArray::New();
		scalars->SetNumberOfComponents(3);
		scalars->SetNumberOfTuples(s->MP_vis.size());
		for (int i=0; i<s->MP_vis.size(); i++)
		{
			scalars->InsertTuple3(i,s->MP_vis(i,0)*255, 0, s->MP_vis_2(i,0)*255);
		}
		vis_polydata->GetPointData()->SetScalars(scalars);//scalars->Delete();//uncomment if there are curvature values associated with polydata

}
void SimGast::shell_cut2polydata1(const shp_shell *s, vtkPolyData* polydata_in,vtkPolyData* polydata_out,vtkPolyData* polydata_fill)
{
	// configure the polydata_in object for the inner shell
	vtkCellArray *_polys1 = vtkCellArray::New();
	for (int i=0; i<this->s->Fc.size(); i++)
	{
		_polys1->InsertNextCell(3);
		_polys1->InsertCellPoint(this->s->Fc[i][0]);
		_polys1->InsertCellPoint(this->s->Fc[i][1]);
		_polys1->InsertCellPoint(this->s->Fc[i][2]);
	}
	vtkPoints* points1 = vtkPoints::New();
	points1->SetDataTypeToDouble();
	points1->SetNumberOfPoints(s->Xin.size());
	for (int i=0; i<s->Xin.size(); i++){	points1->SetPoint(i, s->Xin[i][0], s->Xin[i][1], s->Xin[i][2]);}
	polydata_in->SetPoints(points1);points1->Delete();// Assign points and cells
	polydata_in->SetPolys(_polys1);//polys1->Delete();
	
	// configure the polydata_out object for the outer shell
	vtkPoints* points2 = vtkPoints::New();
	points2->SetDataTypeToDouble();
	points2->SetNumberOfPoints(s->Xout.size());
	for (int i=0; i<s->Xout.size(); i++){	points2->SetPoint(i, s->Xout[i][0], s->Xout[i][1], s->Xout[i][2]);}
	polydata_out->SetPoints(points2);points2->Delete();// Assign points and cells
	polydata_out->SetPolys(_polys1);//polys1->Delete(); // re-use the same mesh from the inner shell

	// configure the polydata_fill object for the fill-in
	vtkCellArray *_polys_fill = vtkCellArray::New();
	for (int i=0; i<this->s->Fcp.size(); i++)
//	for (int i=0; i<25; i++)
	{
		_polys_fill->InsertNextCell(4);						// not given as a triangular mesh
		_polys_fill->InsertCellPoint(this->s->Fcp[i][0]);
		_polys_fill->InsertCellPoint(this->s->Fcp[i][1]);
		_polys_fill->InsertCellPoint(this->s->Fcp[i][2]);
		_polys_fill->InsertCellPoint(this->s->Fcp[i][3]);
	}
	vtkPoints* points3 = vtkPoints::New();
	points3->SetDataTypeToDouble();
	points3->SetNumberOfPoints(s->Xcut.size());
	for (int i=0; i<s->Xcut.size(); i++){	points3->SetPoint(i, s->Xcut[i][0], s->Xcut[i][1], s->Xcut[i][2]);}
	polydata_fill->SetPoints(points3);points3->Delete();// Assign points and cells
	polydata_fill->SetPolys(_polys_fill);//

	// add the scalar fields if needed
	
	if(this->ui.comboBox_sf1->currentIndex()>0 || this->ui.comboBox_sf2->currentIndex()>0 ||this->ui.comboBox_sf3->currentIndex()>0 )
	{
		_min_sf_1 = 0.0;
		_max_sf_1 = 0.0;
		int nsize = 0;
		if(this->ui.comboBox_sf1->currentIndex()>0)
		{
			float tmp;
			for (int i=0; i<s->sf1_vis.size(); i++)
			{
 				if((this->s->sf1_vis[i])>_max_sf_1){_max_sf_1 = this->s->sf1_vis[i];}
				if((this->s->sf1_vis[i])<_min_sf_1){_min_sf_1 = this->s->sf1_vis[i];}
			}
			nsize = s->sf1_vis.size();
		}
		if(this->ui.comboBox_sf2->currentIndex()>0)
		{
			float tmp;
			for (int i=0; i<s->sf2_vis.size(); i++)
			{
 				if((this->s->sf2_vis[i])>_max_sf_1){_max_sf_1 = this->s->sf2_vis[i];}
				if((this->s->sf2_vis[i])<_min_sf_1){_min_sf_1 = this->s->sf2_vis[i];}
			}
			nsize = s->sf2_vis.size();
		}
		if(this->ui.comboBox_sf3->currentIndex()>0)
		{
			float tmp;
			for (int i=0; i<s->sf3_vis.size(); i++)
			{
 				if((this->s->sf3_vis[i])>_max_sf_1){_max_sf_1 = this->s->sf3_vis[i];}
				if((this->s->sf3_vis[i])<_min_sf_1){_min_sf_1 = this->s->sf3_vis[i];}
			}
			nsize = s->sf3_vis.size();
		}
		vtkUnsignedCharArray *scalars = vtkUnsignedCharArray::New();
		scalars->SetNumberOfComponents(3);
		scalars->SetNumberOfTuples(nsize  );	//assuming all scalar fields to hold the same number of elements
		//double fac = 255/_max_sf_1;
		double fac = 255.0/(_max_sf_1-_min_sf_1);
		//std::ostringstream os;
		//os<<"min and max field values displayed: "<<_min_sf_1<<"\t"<<_max_sf_1<<"\n";
		//os<<"Values entered:\n";
		for (int i=0; i<nsize; i++)
		{
			double s1 = 0.0;double s2 = 0.0; double s3 = 0.0;
			if(this->ui.comboBox_sf1->currentIndex()>0){s1 = (this->s->sf1_vis[i] - _min_sf_1)*fac;}
			if(this->ui.comboBox_sf2->currentIndex()>0){s2 = (this->s->sf2_vis[i] - _min_sf_1)*fac;}
			if(this->ui.comboBox_sf3->currentIndex()>0){s3 = (this->s->sf3_vis[i] - _min_sf_1)*fac;}
			scalars->InsertTuple3(i,s1, s2, s3);
			//if(s1>=254.0){os<<i<<"\t"<<s1<<"\t"<<this->s->sf1_vis[i]<<"\t"<<s3<<"\n";}
		}
		polydata_out->GetPointData()->SetScalars(scalars);//scalars->Delete();//uncomment if there are curvature values associated with polydata
		polydata_in->GetPointData()->SetScalars(scalars);
	}
	/*if(this->ui.comboBox_sf1->currentIndex()>0)
	{
		// add the scalar field 1 
		vtkFloatArray *scalars = vtkFloatArray::New();
		float tmp;
		_min_sf_1 = 0.0;
		_max_sf_1 = 0.0;
		for (int i=0; i<s->sf1_vis.size(); i++){
	 		if(abs(this->s->sf1_vis[i])>_max_sf_1){_max_sf_1 = this->s->sf1_vis[i];}
			if((this->s->sf1_vis[i])<_min_sf_1){_min_sf_1 = this->s->sf1_vis[i];}
		}
		for (int i=0; i<s->sf1_vis.size(); i++){scalars->InsertNextValue(this->s->sf1_vis[i]);}
		polydata_out->GetPointData()->SetScalars(scalars);//scalars->Delete();//uncomment if there are curvature values associated with polydata
		polydata_in->GetPointData()->SetScalars(scalars);
	}
*/

}
void SimGast::update_vtk()
{
	if(current_view_type ==0){this->update_vtk1();}		// show cutplane surface
	if(current_view_type ==1){this->update_vtk2();}		// show whole morphology
	if(current_view_type ==2){this->update_vtk3();}		// show whole morphology with active region
	if(current_view_type ==3){this->update_vtk4();}		// show whole morphology with curvature map
}
void SimGast::update_vtk1()
{
	//gets called when the shell cut-plane view needs to be repainted
	this->s->update_cut_plane();
	vtkPolyData*		polydata_in	= vtkPolyData::New();
	vtkPolyData*		polydata_out= vtkPolyData::New();
	vtkPolyData*		polydata_fill= vtkPolyData::New();
	shell_cut2polydata1(this->s, polydata_in, polydata_out, polydata_fill);
	update_display_fields();
	Mapper1_in->SetInputData(polydata_in);
	Mapper1_in->SetScalarRange(_min_sf_1, _max_sf_1);
	Actor1_in->SetMapper( Mapper1_in);

	Mapper1_out->SetInputData(polydata_out);
	Mapper1_out->SetScalarRange(_min_sf_1, _max_sf_1);
	Actor1_out->SetMapper( Mapper1_out);

	Mapper1_fill->SetInputData(polydata_fill);
	Mapper1_fill->SetScalarRange(_min_sf_1, _max_sf_1);
	Actor1_fill->SetMapper( Mapper1_fill);
	//Actor1_fill->GetProperty()->SetOpacity(0.7);

	// test --- add the edges
	vtkExtractEdges *edges = vtkExtractEdges::New();
	edges->SetInputData(polydata_fill);
	vtkPolyDataMapper *edge_mapper = vtkPolyDataMapper::New();
	edge_mapper->SetInputData(edges->GetOutput());
	vtkActor *edge_actor = vtkActor::New();
	edge_actor->SetMapper(edge_mapper);
	edge_actor->GetProperty()->SetColor(1,0.5,0.5);
///

	_ren1->RemoveAllViewProps();
	if(this->axis_on_1){_ren1->AddActor(this->_axes_actor1);}
	if(this->cube_axis_actor1_on)
	{
	// sosi -- add the cube axes actor
	this->_cube_axes_actor1->SetBounds(Actor1_out->GetBounds());
	this->_cube_axes_actor1->SetCamera(_ren1->GetActiveCamera());
	_ren1->AddActor(this->_cube_axes_actor1);
	}
	if(this->vm_actor_on){_ren1->AddActor(this->vmActor);}
	_ren1->AddActor(this->Actor1_in);
	_ren1->AddActor(this->Actor1_out);
	_ren1->AddActor(this->Actor1_fill);
	_ren1->AddActor(edge_actor);
	_ren1->ResetCamera();
	_ren1->GetRenderWindow()->Render();
}
void SimGast::update_vtk2()
{
	//gets called when the shell needs to be repainted with the whole morphology of the mid-plane
	this->s->update();
	//vtkPolyData*		polydata	= vtkPolyData::New();
	shell2polydata1();
	update_display_fields();
	Mapper1_in->SetInputData(vis_polydata);
	//Mapper1_in->SetScalarRange(_min_sf_1, _max_sf_1);
	Actor1_in->SetMapper( Mapper1_in);

	_ren1->RemoveAllViewProps();
	if(this->axis_on_1){_ren1->AddActor(this->_axes_actor1);}
	if(this->cube_axis_actor1_on)
	{
	// sosi -- add the cube axes actor
	this->_cube_axes_actor1->SetBounds(Actor1_in->GetBounds());
	this->_cube_axes_actor1->SetCamera(_ren1->GetActiveCamera());
	_ren1->AddActor(this->_cube_axes_actor1);
	}

	if(this->vm_actor_on){_ren1->AddActor(this->vmActor);}
	_ren1->AddActor(this->Actor1_in);
	_ren1->ResetCamera();
	_ren1->GetRenderWindow()->Render();
}

void SimGast::update_vtk3()
{
	//show whole morphology with active region
	this->s->update_scalar_fields();
	this->s->update();
	shell2polydata1_active_region();
	update_display_fields();
	Mapper1_in->SetInputData(vis_polydata);
	//Mapper1_in->SetScalarRange(_min_sf_1, _max_sf_1);
	Actor1_in->SetMapper( Mapper1_in);

	_ren1->RemoveAllViewProps();
	if(this->axis_on_1){_ren1->AddActor(this->_axes_actor1);}
	if(this->cube_axis_actor1_on)
	{
		// sosi -- add the cube axes actor
		this->_cube_axes_actor1->SetBounds(Actor1_in->GetBounds());
		this->_cube_axes_actor1->SetCamera(_ren1->GetActiveCamera());
		_ren1->AddActor(this->_cube_axes_actor1);
	}
	if(this->vm_actor_on){_ren1->AddActor(this->vmActor);}
	_ren1->AddActor(this->Actor1_in);
	_ren1->ResetCamera();
	_ren1->GetRenderWindow()->Render();
}
void SimGast::update_vtk4()
{
	_ren1->RemoveAllViewProps();
	//show whole morphology with active region
	this->s->update();
	this->smv->surfaceGenX(this->s->srf_m->xc, this->s->srf_m->yc, this->s->srf_m->zc); // use fast surface generation for X only
	//std::ostringstream os;os<<"--------------------"<<this->smv->Xu<<std::endl;report(os.str());os.str("");
	vtkPoints* points = vtkPoints::New();
	points->SetDataTypeToDouble();
	points->SetNumberOfPoints(this->smv->n_points);
	for (int i=0; i<this->smv->n_points; i++)
	{	
		points->SetPoint(i, this->smv->Xu(i,0), this->smv->Xu(i,1), this->smv->Xu(i,2));/// fill in the points
	}
	vis_polydata->SetPoints(points);points->Delete();// Assign points
	////////////////////////////////////////////////


	///////////////// plot the local mean curvature ///////////////
	this->smv->_curv_calc = 1;
	this->smv->update_tri();
	_min_sf_1 = 0.0;
	_max_sf_1 = 0.0;
	vtkFloatArray *scalars = vtkFloatArray::New();
	for (int i=0; i<smv->n_points; i++)	
	{
		if((this->smv->H[i])>_max_sf_1){_max_sf_1 = this->smv->H[i];}
		if((this->smv->H[i])<_min_sf_1){_min_sf_1 = this->smv->H[i];}
		scalars->InsertNextValue(this->smv->H[i]);
	}
	vis_polydata->GetPointData()->SetScalars(scalars);//scalars->Delete();
	/////////////////////////////////////////
	optimMapper->SetInputData(vis_polydata);
	optimMapper->SetScalarRange(_min_sf_1, _max_sf_1);
	optimActor->SetMapper( optimMapper);

	////////////// add the scalar bar
	scalarBar->SetLookupTable(optimMapper->GetLookupTable());
	_ren1->AddActor(scalarBar);
	///////////////
	if(this->axis_on_1){_ren1->AddActor(this->_axes_actor1);}
	if(this->cube_axis_actor1_on)
	{
		// sosi -- add the cube axes actor
		this->_cube_axes_actor1->SetBounds(Actor1_in->GetBounds());
		this->_cube_axes_actor1->SetCamera(_ren1->GetActiveCamera());
		_ren1->AddActor(this->_cube_axes_actor1);
	}
	if(this->vm_actor_on){_ren1->AddActor(this->vmActor);}
	_ren1->AddActor(this->optimActor);
	_ren1->ResetCamera();
	_ren1->GetRenderWindow()->Render();
	this->ui.vtkWidget1->update();
}

void SimGast::update_display_fields()
{
	////////////////
	ui.lineEdit_Lmax->setText(QString::number(double(this->Lmax)));
	ui.lineEdit_gdim->setText(QString::number(double(this->gdim)));
	ui.lineEdit_D->setText(QString::number(double(s->get_D())));
	ui.lineEdit_VMscale->setText(QString::number(double(s->get_VMscale())));
	ui.lineEdit_gLmax->setText(QString::number(double(this->max_gL_max)));
	ui.lineEdit_young->setText(QString::number(double(s->get_Young())));
	ui.lineEdit_poisson->setText(QString::number(double(s->get_Poisson())));


	ui.lineEdit_gamma1->setText(QString::number(double(s->get_psnl())));
	ui.lineEdit_gamma2->setText(QString::number(double(s->get_ptw())));
	ui.lineEdit_gamma3->setText(QString::number(double(s->get_phkb())));
	ui.lineEdit_mp_cutoff->setText(QString::number(double(s->get_mp_cutoff())));
	ui.lineEdit_stiff_fac->setText(QString::number(double(s->get_mp_fac())));
	ui.lineEdit_active_region_curvature->setText(QString::number(double(s->get_mp_co())));

	ui.lineEdit_gamma1_2->setText(QString::number(double(s->get_psf1_2())));
	ui.lineEdit_gamma2_2->setText(QString::number(double(s->get_psf2_2())));
	ui.lineEdit_gamma3_2->setText(QString::number(double(s->get_psf3_2())));
	ui.lineEdit_mp_cutoff_2->setText(QString::number(double(s->get_mp_cutoff_2())));
	ui.lineEdit_stiff_fac_2->setText(QString::number(double(s->get_mp_fac_2())));
	ui.lineEdit_active_region_curvature_2->setText(QString::number(double(s->get_mp_co_2())));

	////////////////
}

//////////////////////// FITTING HELPER FUNCTIONS ////////////////////
void SimGast::on_pushButton_start_clicked()
{
	stop_optimization=false;
	if(this->ui.comboBox_sf1->currentIndex()>0 & this->ui.comboBox_sf2->currentIndex()>0 & this->ui.comboBox_sf3->currentIndex()>0)
	{
		// prepare the shell morphology for the limited Lmax specified for fitting
		const QString text1  = this->ui.lineEdit_PSO_Lmax_fit->text();
		maxLfit = text1.toInt();
		this->s->srf_m->flush_after_L(maxLfit);
		this->s->srf_u->flush_after_L(maxLfit);
		//
		std::ostringstream os;os<<"---------STARTING CONFIGURATION --------------"<<*s;//report(os.str());
		/*// sosi-- modify a coefficient and calculate energy
		this->s->srf_m->update();
		
		os<<"Value before\n"<<this->s->srf_m->xc(1,0)<<std::endl;
		this->s->srf_m->xc(1,0) = this->s->srf_m->xc(1,0) *1.5;
		this->s->srf_m->update();
		os<<"Value after\n"<<this->s->srf_m->xc(1,0)<<"\n"<<*s<<std::endl;
		//this->s->srf_m->update();
		*/
		int old_tri_n = this->s->srf_m->sm->tri_n;
		this->s->srf_m->set_new_spherical_mesh(tri_n_intersection_test);		// we use a low resolution mesh for intersection tests
		double E = this->s->get_energy(1);
		os<<"\n"<<"Active region area fraction:"<<s->active_area_fraction<<std::endl;
		//os<<"Vout Vin Aout Ain:"<<<<"\n"<<s->vin<<"\t"<<s->aout<<"\t"<<s->ain<<"\t"<<std::endl;
		os<<"Volume enclosed by outer surface\t"<<s->vout<<std::endl;
		os<<"Area of outer surface\t"<<s->aout<<std::endl;
		os<<"Volume enclosed by inner surface\t\t"<<s->vin<<std::endl;
		os<<"Area of inner surface\t"<<s->ain<<std::endl;

		os<<"Using input morphology both as reference and starting configurations."<<std::endl;
		os<<"Shear and stretch: "<<this->s->E_nHook<<"\tBending: "<<this->s->Eb<<"\nTotal energy:\n"<<E<<std::endl;
		report(os.str());os.str("");

		this->mc_start();
		if(stop_optimization){os<<"Optimization stopped at user's request";report(os.str());os.str("");}
		this->s->srf_m->set_new_spherical_mesh(old_tri_n);		// reset the spherical mesh to what it was before
	}
	else
	{
		this->message("Please select three scalar fields.");
	}

}
void SimGast::on_pushButton_start_PSO_clicked()
{
	stop_optimization=false;
	//if(this->ui.comboBox_sf1->currentIndex()>0 & this->ui.comboBox_sf2->currentIndex()>0 & this->ui.comboBox_sf3->currentIndex()>0)
	//{
		// prepare the shell morphology for the limited Lmax specified for fitting
		const QString text1  = this->ui.lineEdit_PSO_Lmax_fit->text();
		maxLfit = text1.toInt();
		this->s->srf_m->flush_after_L(maxLfit);
		this->s->srf_u->flush_after_L(maxLfit);
		//
		std::ostringstream os;os<<"---------STARTING CONFIGURATION --------------"<<*s;//report(os.str());
		/*// sosi-- modify a coefficient and calculate energy
		this->s->srf_m->update();
		
		os<<"Value before\n"<<this->s->srf_m->xc(1,0)<<std::endl;
		this->s->srf_m->xc(1,0) = this->s->srf_m->xc(1,0) *1.5;
		this->s->srf_m->update();
		os<<"Value after\n"<<this->s->srf_m->xc(1,0)<<"\n"<<*s<<std::endl;
		//this->s->srf_m->update();
		*/
		int old_tri_n = this->s->srf_m->sm->tri_n;
		this->s->srf_m->set_new_spherical_mesh(tri_n_intersection_test);		// we use a low resolution mesh for intersection tests
		double E = this->s->get_energy(1);
		os<<"\n"<<"Active region area fraction:"<<s->active_area_fraction<<std::endl;
		//os<<"Vout Vin Aout Ain:"<<<<"\n"<<s->vin<<"\t"<<s->aout<<"\t"<<s->ain<<"\t"<<std::endl;
		os<<"Volume enclosed by outer surface\t"<<s->vout<<std::endl;
		os<<"Area of outer surface\t"<<s->aout<<std::endl;
		os<<"Volume enclosed by inner surface\t\t"<<s->vin<<std::endl;
		os<<"Area of inner surface\t"<<s->ain<<std::endl;

		os<<"Using input morphology both as reference and starting configurations."<<std::endl;
		os<<"Shear and stretch: "<<this->s->E_nHook<<"\tBending: "<<this->s->Eb<<"\nTotal energy:\n"<<E<<std::endl;
		report(os.str());os.str("");

		this->pso_start();
		if(stop_optimization){os<<"Optimization stopped at user's request";report(os.str());os.str("");}
		this->s->srf_m->set_new_spherical_mesh(old_tri_n);		// reset the spherical mesh to what it was before
	//}
	//else
	//{
	//	this->message("Please select three scalar fields.");
	//}

}
void SimGast::generate_vm_actor_hi_res()
{//// generate the vitteline membrane actor
	this->s->prepare_energy_calc();
	this->smv->surfaceGen(this->s->srf_u->xc, this->s->srf_u->yc, this->s->srf_u->zc); // 
	this->smv->update_tri();
	this->smv->sfGen(this->s->Dc);	// generates the scalar field corresponding to the thickness
	// generate the inner and outer planes
	double fac = 10000;
	vector < vector<double> > Xvm;
	int n_points = this->smv->n_points;
	//Xin.resize(n_points);
	Xvm.resize(n_points);
	for(int i=0;i<n_points;i++)
	{
		Xvm[i].resize(3);
		Xvm[i][0] = smv->X[i][0]-this->s->fix_normals_sign*smv->vN[i][0]*(smv->sf[i]/2+s->get_VMscale());
		Xvm[i][1] = smv->X[i][1]-this->s->fix_normals_sign*smv->vN[i][1]*(smv->sf[i]/2+s->get_VMscale());
		Xvm[i][2] = smv->X[i][2]-this->s->fix_normals_sign*smv->vN[i][2]*(smv->sf[i]/2+s->get_VMscale());
	}

	vtkPoints* vmpoints = vtkPoints::New();
	vmpoints->SetDataTypeToDouble();
	vmpoints->SetNumberOfPoints(this->smv->n_points);
	for (int i=0; i<this->smv->n_points; i++)
	{	
		vmpoints->SetPoint(i,Xvm[i][0],Xvm[i][1], Xvm[i][2]);/// fill in the points
	}
	vm_polydata->SetPoints(vmpoints);vmpoints->Delete();// Assign points
	vtkDepthSortPolyData* ds = vtkDepthSortPolyData::New();
	ds->SetInputData(vm_polydata);
	ds->SetDirectionToBackToFront();
	ds->SetVector(1, 1, 1);
	ds->SetCamera(_ren1->GetActiveCamera());
	ds->SortScalarsOn();
	ds->Update();

	vmMapper->SetInputData(ds->GetOutput());
	vmActor->SetMapper( vmMapper);
	vmActor->GetProperty()->SetOpacity(.35);
}
void SimGast::generate_vm_actor()
{
	//this->s->update();
	this->s->prepare_energy_calc();
	//// generate the vitteline membrane actor as in srf_m
	
	vtkCellArray * polys = vtkCellArray::New();
	for (int i=0; i<this->s->srf_m->sm->n_faces; i++)
	{
		polys->InsertNextCell(3);
		polys->InsertCellPoint(this->s->srf_m->sm->f0[i]-1);
		polys->InsertCellPoint(this->s->srf_m->sm->f1[i]-1);
		polys->InsertCellPoint(this->s->srf_m->sm->f2[i]-1);
	}
	vtkPolyData *vm_polydata_vis	= vtkPolyData::New();
	vm_polydata_vis->SetPolys(polys);//polys->Delete();
	/////

	vtkPoints* vmpoints = vtkPoints::New();
	vmpoints->SetDataTypeToDouble();
	vmpoints->SetNumberOfPoints(this->s->srf_m->sm->n_points);
	for (int i=0; i<this->s->srf_m->sm->n_points; i++)
	{	
		vmpoints->SetPoint(i, s->Xvm[i][0], s->Xvm[i][1], s->Xvm[i][2]);/// fill in the points
	}
	vm_polydata_vis->SetPoints(vmpoints);vmpoints->Delete();// Assign points
	vtkDepthSortPolyData* ds = vtkDepthSortPolyData::New();
	ds->SetInputData(vm_polydata_vis);
	ds->SetDirectionToBackToFront();
	ds->SetVector(1, 1, 1);
	ds->SetCamera(_ren1->GetActiveCamera());
	ds->SortScalarsOn();
	ds->Update();

	vmMapper->SetInputData(ds->GetOutput());
	vmActor->SetMapper( vmMapper);
	vmActor->GetProperty()->SetOpacity(.35);

}
void SimGast::update_vis_optim_hi_res()
{	// designed to be fast -- at the moment shows only the mid-plane surface
	//gets called when the structure needs to be repainted during optimization
	// requires prior definition of vis_polys and smv, and updating of smv using smv->surfaceGen(xc, yc, zc);
	// definition of vis_polys is done at initialization
	//this->s->update();	// we assume that update has happened during optimization already, so don't repeat it here


	/*
	// uncomment to use Fast -- just the mid-plane surface
	this->smv->surfaceGenX(this->s->srf_m->xc, this->s->srf_m->yc, this->s->srf_m->zc); // use fast surface generation for X only
	//std::ostringstream os;os<<"--------------------"<<this->smv->Xu<<std::endl;report(os.str());os.str("");
	vtkPoints* points = vtkPoints::New();
	points->SetDataTypeToDouble();
	points->SetNumberOfPoints(this->smv->n_points);
	for (int i=0; i<this->smv->n_points; i++)
	{	
		points->SetPoint(i, this->smv->Xu(i,0), this->smv->Xu(i,1), this->smv->Xu(i,2));/// fill in the points
	}
	vis_polydata->SetPoints(points);points->Delete();
	*/
	/// slow ---- plot the outer surface
	this->smv->surfaceGen(this->s->srf_m->xc, this->s->srf_m->yc, this->s->srf_m->zc); // 
	this->smv->update_tri();
	this->smv->sfGen(this->s->Dc);	// generates the scalar field corresponding to the thickness
	// generate the inner and outer planes
	double fac = 10000;
	/*vector<vector<double>> Xin, Xout;
	int n_points = this->smv->n_points;
	Xin.resize(n_points);
	Xout.resize(n_points);
	for(int i=0;i<n_points;i++)
	{
		Xout[i].resize(3);
		Xout[i][0] = smv->X[i][0]-smv->vN[i][0]*smv->sf[i]/2;
		Xout[i][1] = smv->X[i][1]-smv->vN[i][1]*smv->sf[i]/2;
		Xout[i][2] = smv->X[i][2]-smv->vN[i][2]*smv->sf[i]/2;
		//Xin[i].resize(3);
		//Xin[i][0] = smv->X[i][0]+smv->vN[i][0]*smv->sf[i]/2;
		//Xin[i][1] = smv->X[i][1]+smv->vN[i][1]*smv->sf[i]/2;
		//Xin[i][2] = smv->X[i][2]+smv->vN[i][2]*smv->sf[i]/2;
	}*/
	vtkPoints* points = vtkPoints::New();
	points->SetDataTypeToDouble();
	points->SetNumberOfPoints(this->smv->n_points);
	for (int i=0; i<this->smv->n_points; i++)
	{	
		points->SetPoint(i, smv->X[i][0]-smv->vN[i][0]*smv->sf[i]/2, smv->X[i][1]-smv->vN[i][1]*smv->sf[i]/2, smv->X[i][2]-smv->vN[i][2]*smv->sf[i]/2);/// fill in the points
	}
	vis_polydata->SetPoints(points);points->Delete();
	///////////////// plot the local mean curvature ///////////////
	this->smv->_curv_calc = 1;
	this->smv->update_tri();
	_min_sf_1 = 0.0;
	_max_sf_1 = 0.0;
	vtkFloatArray *scalars = vtkFloatArray::New();
	for (int i=0; i<smv->n_points; i++)	
	{
		if((this->smv->H[i])>_max_sf_1){_max_sf_1 = this->smv->H[i];}
		if((this->smv->H[i])<_min_sf_1){_min_sf_1 = this->smv->H[i];}
		scalars->InsertNextValue(this->smv->H[i]);
	}
	vis_polydata->GetPointData()->SetScalars(scalars);scalars->Delete();
	/////////////////////////////////////////
	optimMapper->SetInputData(vis_polydata);
	optimMapper->SetScalarRange(_min_sf_1, _max_sf_1);
	optimActor->SetMapper( optimMapper);
	/**/
	//////////////////////////////////////////

	_ren1->RemoveAllViewProps();
	
	////////////// add the scalar bar
	scalarBar->SetLookupTable(optimMapper->GetLookupTable());
	scalarBar->SetVisibility(true);
	_ren1->AddActor(scalarBar);
	_ren1->AddActor(this->optimActor);
	if(this->vm_actor_on){_ren1->AddActor(this->vmActor);}
	
	//_ren1->ResetCamera();
	_ren1->GetRenderWindow()->Render();
	
	// cleanup
}

void SimGast::update_vis_optim()
{	
	/// 
	vtkPoints* points = vtkPoints::New();
	points->SetDataTypeToDouble();
	points->SetNumberOfPoints(this->s->srf_m->sm->n_points);
	for (int i=0; i<this->s->srf_m->sm->n_points; i++)
	{	
		points->SetPoint(i, s->Xout[i][0], s->Xout[i][1], s->Xout[i][2]);/// fill in the points
	}
	vtkCellArray * polys = vtkCellArray::New();
	for (int i=0; i<this->s->srf_m->sm->n_faces; i++)
	{
		polys->InsertNextCell(3);
		polys->InsertCellPoint(this->s->srf_m->sm->f0[i]-1);
		polys->InsertCellPoint(this->s->srf_m->sm->f1[i]-1);
		polys->InsertCellPoint(this->s->srf_m->sm->f2[i]-1);
	}
	vtkPolyData *vm_polydata_vis	= vtkPolyData::New();
	vm_polydata_vis->SetPolys(polys);polys->Delete();
	vm_polydata_vis->SetPoints(points);points->Delete();




	///////////////// plot the local mean curvature ///////////////
	this->s->srf_m->sm->_curv_calc = 1;
	this->s->srf_m->sm->update_tri();
	_min_sf_1 = 0.0;
	_max_sf_1 = 0.0;
	vtkFloatArray *scalars = vtkFloatArray::New();
	for (int i=0; i<s->srf_m->sm->n_points; i++)	
	{
		if((this->s->srf_m->sm->H[i])>_max_sf_1){_max_sf_1 = this->s->srf_m->sm->H[i];}
		if((this->s->srf_m->sm->H[i])<_min_sf_1){_min_sf_1 = this->s->srf_m->sm->H[i];}
		scalars->InsertNextValue(this->s->srf_m->sm->H[i]);
	}
	vm_polydata_vis->GetPointData()->SetScalars(scalars);scalars->Delete();
	/////////////////////////////////////////
	optimMapper->SetInputData(vm_polydata_vis);
	optimMapper->SetScalarRange(_min_sf_1, _max_sf_1);
	optimActor->SetMapper( optimMapper);
	/**/
	//////////////////////////////////////////

	_ren1->RemoveAllViewProps();
	
	////////////// add the scalar bar
	scalarBar->SetLookupTable(optimMapper->GetLookupTable());
	scalarBar->SetVisibility(true);
	_ren1->AddActor(scalarBar);
	_ren1->AddActor(this->optimActor);
	if(this->vm_actor_on){_ren1->AddActor(this->vmActor);}
	//_ren1->ResetCamera();
	_ren1->GetRenderWindow()->Render();
}


void SimGast::on_pushButton_stop_clicked()
{
	stop_optimization = true;

}
void SimGast::on_pushButton_stop_PSO_clicked()
{
	stop_optimization = true;

}
void SimGast::on_pushButton_stop_COBYLA_clicked()
{
	stop_optimization = true;

}
void SimGast::on_pushButton_stop_SBPLX_clicked()
{
	stop_optimization = true;

}
void SimGast::on_pushButton_stop_ISRES_clicked()
{
	stop_optimization = true;

}
int SimGast::intersection_tests()
	// call only after prepare_energy_calc() has been called at some point to define the VM
	// and after update() has been called to update all surfaces
{
	int intersection = 1;
	if(s->VM_intersect()==0)
	{
		//intersection = 3;
		//if(s->self_intersect()==0)
		//{
		//	intersection = 0;
		//}
		if(s->srf_m->self_intersect()==0)
		{
			intersection = 0;
		}
	}
	return intersection;
}

int SimGast::vtk_intersection_tests()
	// call only after prepare_energy_calc() has been called at some point to define the VM
	// and after update() has been called to update all surfaces
{
	int intersection = 0;
	int coplanar;
    double surfaceid[2];
    double tolerance = 0.00001;
    double isectpt1[3], isectpt2[3];
	double* V0;
	double* V1;
	double* V2;
	double* U0;
	double* U1;
	double* U2;

	if(this->opt_intersection_test_vm)
	{
		for(int i = 0;i<s->srf_m->sm->n_faces;i++)	// test for intersection with vitteline membrane
		{
			V0 = (&s->Xout[int(s->srf_m->sm->f0[i]-1)][0]);
			V1 = (&s->Xout[int(s->srf_m->sm->f1[i]-1)][0]);
			V2 = (&s->Xout[int(s->srf_m->sm->f2[i]-1)][0]);
			for(int j = 0;j<s->srf_m->sm->n_faces;j++)
			{
				U0 = (&s->Xvm[int(s->srf_m->sm->f0[j]-1)][0]);
				U1 = (&s->Xvm[int(s->srf_m->sm->f1[j]-1)][0]);
				U2 = (&s->Xvm[int(s->srf_m->sm->f2[j]-1)][0]);
				//intersection = vtkIntersectionPolyDataFilter::TriangleTriangleIntersection(V0,V1,V2,U0,U1,U2, coplanar,isectpt1,isectpt2, surfaceid, tolerance);
				if(intersection){return 1;}
			}
		}
	}
	if(this->opt_intersection_test_inner)
	{
		if(intersection==0)		// then test for self-intersection of inner surface of the shell
		{
			for(int i = 0;i<s->srf_m->sm->n_faces;i++)
			{
				V0 = (&s->Xin[int(s->srf_m->sm->f0[i]-1)][0]);
				V1 = (&s->Xin[int(s->srf_m->sm->f1[i]-1)][0]);
				V2 = (&s->Xin[int(s->srf_m->sm->f2[i]-1)][0]);
				for(int j = i+1;j<s->srf_m->sm->n_faces;j++)
				{
					U0 = (&s->Xin[int(s->srf_m->sm->f0[j]-1)][0]);
					U1 = (&s->Xin[int(s->srf_m->sm->f1[j]-1)][0]);
					U2 = (&s->Xin[int(s->srf_m->sm->f2[j]-1)][0]);
					if(!share_vert((s->srf_m->sm->f0[i]-1),(s->srf_m->sm->f1[i]-1),(s->srf_m->sm->f2[i]-1),(s->srf_m->sm->f0[j]-1),(s->srf_m->sm->f1[j]-1),(s->srf_m->sm->f2[j]-1)))
					{	
						//intersection = vtkIntersectionPolyDataFilter::TriangleTriangleIntersection(V0,V1,V2,U0,U1,U2, coplanar,isectpt1,isectpt2, surfaceid, tolerance);
					}
					if(intersection){return 1;}
				}
			}
		}
	}
	if(this->opt_intersection_test_outer)
	{
		if(intersection==0)		// then test for self-intersection of outer surface of the shell
		{
			for(int i = 0;i<s->srf_m->sm->n_faces;i++)
			{
				V0 = (&s->Xout[int(s->srf_m->sm->f0[i]-1)][0]);
				V1 = (&s->Xout[int(s->srf_m->sm->f1[i]-1)][0]);
				V2 = (&s->Xout[int(s->srf_m->sm->f2[i]-1)][0]);
				for(int j = i+1;j<s->srf_m->sm->n_faces;j++)
				{
					U0 = (&s->Xout[int(s->srf_m->sm->f0[j]-1)][0]);
					U1 = (&s->Xout[int(s->srf_m->sm->f1[j]-1)][0]);
					U2 = (&s->Xout[int(s->srf_m->sm->f2[j]-1)][0]);
					if(!share_vert((s->srf_m->sm->f0[i]-1),(s->srf_m->sm->f1[i]-1),(s->srf_m->sm->f2[i]-1),(s->srf_m->sm->f0[j]-1),(s->srf_m->sm->f1[j]-1),(s->srf_m->sm->f2[j]-1)))
					{	
						//intersection = vtkIntersectionPolyDataFilter::TriangleTriangleIntersection(V0,V1,V2,U0,U1,U2, coplanar,isectpt1,isectpt2, surfaceid, tolerance);
					}
					if(intersection){return 1;}
				}
			}
		}
	}
	return intersection;
}

void SimGast::generate_XYplot(std::vector<double> y)
{
    
	if(y.size()>1)
	{
	int numPoints = y.size();
			vtkTable* table = vtkTable::New();
		table->SetNumberOfRows(numPoints);
		vtkFloatArray * arrX = vtkFloatArray::New();
		vtkFloatArray * arrY = vtkFloatArray::New();
		arrX->SetName("iteration number");
		arrY->SetName("Energy");
		arrX->SetNumberOfValues(numPoints);
		arrY->SetNumberOfValues(numPoints);
		for(int i = 0;i<numPoints;i++)
		{
			arrX->SetValue(i,i);
			arrY->SetValue(i,y[i]);
		}
		table->AddColumn(arrX);
		table->AddColumn(arrY);
		vtkChartXY* chart = vtkChartXY::New();
		XYview->GetScene()->ClearItems();
		XYview->GetScene()->AddItem(chart);
		chart->ClearPlots();
		chart->GetAxis(vtkAxis::LEFT)->SetTitle("Energy");
		chart->GetAxis(vtkAxis::BOTTOM)->SetTitle("iteration number");
		vtkPlot* curve = chart->AddPlot(vtkChart::LINE);
		curve->SetInputData(table, 0, 1);
		curve->SetColor(0, 255, 0, 255);
		curve->SetWidth(1.0);
		XYview->Render();
		// cleanup
		table->Delete();
		arrX->Delete();
		arrY->Delete();
		chart->Delete();

	}
     /**/
}

//////////////////////// FITTING --- MMC////////////////////
void SimGast::change_configuration(int &coord, int &cix, double &old_val)
{
	//this->s->srf_m->xc(1,0) = this->s->srf_m->xc(1,0) *1.5;
	int maxix = (this->maxLfit + 1)*(this->maxLfit + 1);
	//int maxix = this->s->srf_m->xc.size();
	coord = int(unifRand(3));// generate a random integer between 1 and 3 to decide on the coordinate that will be changed x y or z
	cix = int(unifRand(maxix-1)); // we are avoiding change of the first coefficient
	double new_val = 0.0;

	if(coord==1)
	{
		old_val = this->s->srf_m->xc(cix,0);
		new_val = randn_notrig(this->s->srf_m->xc(cix,0), this->MCsig[cix]);
		this->s->srf_m->xc(cix,0) = new_val;
	}
	if(coord==2)
	{
		old_val = this->s->srf_m->yc(cix,0);
		new_val = randn_notrig(this->s->srf_m->yc(cix,0), this->MCsig[cix]);
		this->s->srf_m->yc(cix,0) = new_val;
	}
	if(coord==3)
	{
		old_val = this->s->srf_m->zc(cix,0);
		new_val = randn_notrig(this->s->srf_m->zc(cix,0), this->MCsig[cix]);
		this->s->srf_m->zc(cix,0) = new_val;
	}

	
	//std::ostringstream os;os<<"\tcoord: "<<coord<<"\t"<<"cix: "<<cix<<"\nOld val\tNew val\n"<<old_val<<"\t"<<new_val<<std::endl;report(os.str());

}
void SimGast::mc_start()
{
	this->generate_vm_actor_hi_res();
	_ren1->ResetCamera();
	std::ostringstream os;os<<"---------STARTING METROPOLIS MONTECARLO OPTIMIZATION-----------"<<std::endl;
	const QString text1  = this->ui.lineEdit_MMC_Lmax_fit->text();
	maxLfit = text1.toInt();
	const QString text2  = this->ui.lineEdit_MMC_maxiter->text();
	int niter = text2.toInt();
	const QString text3  = this->ui.lineEdit_MMC_sigma->text();
	double sig = text3.toDouble();
	const QString text4  = this->ui.lineEdit_MMC_TolCon->text();
	//double TolCon = this->s->srf_m->V *text4.toDouble()/100.0;
	double volPen = text4.toDouble();
	const QString text5  = this->ui.lineEdit_MMC_T->text();
	double fac1 = 1.0/text5.toDouble();// temperature factor
    //double volPen = 1e3;		// Penalty factor for volume constraint violation
	//////////
	os<<"Maximum number of iterations: "<<niter<<"\nConstraint penalty factor: "<<volPen<<"\n";
	os<<"\nFitting coefficients up to Lmax: "<<maxLfit<<"\n"<<"\nTemperature: "<<1.0/fac1<<"\n"<<"------------------------------";report(os.str());os.str("");
	this->MCsig.resize(this->s->srf_m->xc.size());
	for(int i = 0;i<this->s->srf_m->xc.size();i++){this->MCsig[i] = sig;}	// fill MCsig with sig values
	int intersection = this->vtk_intersection_tests();
	os<<"Starting-shape intersection test result: "<<intersection<<std::endl;report(os.str());os.str("");
	std::vector<double> Energy_vec;
	Energy_vec.reserve(niter);
	if(intersection==0)		// if the starting shape is good then we will use it
	{
		
		double Eo = this->s->get_energy(1);	// this is the energy at the beginning (currently accepted)
		double Vo = this->s->srf_u->V;		// set the volume constraint to that of the starting configuration
		//MatrixXd xco = this->s->srf_m->xc;	// currently accepted configuration
		//MatrixXd yco = this->s->srf_m->yc;
		//MatrixXd zco = this->s->srf_m->zc;
		/// for monitoring the MMC:
		int nacc = 0;		// total number of acceptances 
		//int nacc_prob = 0;	// total number of acceptances after probability test
		int nrej_inter = 0;		// total number of rejections due to intersection
		int nrej_no_inter = 0;	// total number of rejections excluding intersection rejections
		int feval_count = 0;  // total count of function evaluations (within a monitoring cycle)
		int mon_cyc  = 100;  // monitoring cycle steps

		////////////////////////////////////
		double E = 0.0;
		double old_val = 0.0;		// here we will store the coefficient value before change
		int coord = 0;				// which coordinate was changed x,y,z
		int cix = 0;				// which coefficient was changed (linear index)
		
		///// the actual MMC iterations start here
		/////////////////////////////////////////
		stop_optimization = false;
		int iter = 0;
		bool volCon = 0;
		double Evol = 0.0;
		this->update_vis_optim_hi_res();		// let's look at it before the iterations start
		Energy_vec.push_back(Eo);
		while(iter<niter & stop_optimization==0)
		{
			QApplication::processEvents(QEventLoop::AllEvents);		//sosi -- in the future start a separate thread for the optimization
			//os<<"iter:\t"<<iter<<"\t";
			//os<<std::endl<<*s->srf_m<<std::endl;report(os.str());os.str("");
			this->change_configuration(coord, cix, old_val);					// randomly change the configuration
			//this->s->srf_m->update();							   // update the geometry and surface derivatives -- this is all we need at the moment
			this->s->update();
			intersection = this->vtk_intersection_tests();
			if(intersection==0)								// check intersection constraint
			{
				//os<<"\tintersection:\t"<<intersection<<"\tdV:\t"<<std::abs(Vo-this->s->srf_m->V)<<std::endl;
				//volCon = std::abs(Vo-this->s->srf_m->V)<TolCon;
				//if(volCon)	// check the volume constraint
				//{
					E = this->s->get_energy();
					Evol = volPen * std::abs(Vo-this->s->srf_m->V);
					E = E + Evol;
					feval_count++;
					//os<<"iterinner :\t"<<iter<<"\t"<<"Eold: "<<Eo<<"\tEtry: "<<E<<std::endl;
					//this->update_vis_optim();
					//this->update_vtk1();
					//update_vtk();
					if(Eo<E)			// the case when we might change the configuration depending on how much bigger E is.
					{
						double p = std::exp(fac1*(Eo-E));
						if(p<double(rand()/double(RAND_MAX)))
						{
							//this->update_vis_optim();
							if(coord==1){this->s->srf_m->xc(cix,0) = old_val;}
							if(coord==2){this->s->srf_m->yc(cix,0) = old_val;}
							if(coord==3){this->s->srf_m->zc(cix,0) = old_val;}
							nrej_no_inter++;
						}
						else
						{
							Eo = E;
							nacc++;
							Energy_vec.push_back(Eo);
							this->update_vis_optim_hi_res();//os<<" --- probability accept"<<fac1<<"\tProbability: "<<p;
							this->generate_XYplot(Energy_vec);
						}
					}
					if(Eo>E)			// the case when we will automatically accept
					{
						Eo = E;
						nacc++;
						Energy_vec.push_back(Eo);
						this->update_vis_optim_hi_res();//os<<" --- accept";
						this->generate_XYplot(Energy_vec);
					}
				//}
			}else
			{
				//this->update_vis_optim();

				if(coord==1){this->s->srf_m->xc(cix,0) = old_val;}
				if(coord==2){this->s->srf_m->yc(cix,0) = old_val;}
				if(coord==3){this->s->srf_m->zc(cix,0) = old_val;}
				nrej_inter++;
				os<<"intersection :\t"<<intersection;report(os.str());os.str("");
			}
			os<<"Iteration\t"<<"Intersect\t"<<"Volcon %\t"<<"Accept\t"<<"Reject\t"<<"FunEvals\n";
			os<<iter+1<<"\t"<<intersection<<"\t"<<this->s->srf_m->V/Vo*100<<"\t"<<nacc<<"\t"<<(nrej_no_inter + nrej_inter)<<"\t"<<feval_count;//report(os.str());os.str("");
			report(os.str());os.str("");
			iter++;
		}
	}
}
//////////////////////// FITTING --- PSO ////////////////
void SimGast::pso_start()
{
	this->generate_vm_actor_hi_res();
	_ren1->ResetCamera();
	std::ostringstream os;os<<"---------STARTING PARTICLE SWARM OPTIMIZATION-----------"<<std::endl;
	const QString text1  = this->ui.lineEdit_PSO_Lmax_fit->text();
	maxLfit = text1.toInt();
	const QString text2  = this->ui.lineEdit_PSO_maxiter->text();
	int niter = text2.toInt();
	const QString text3  = this->ui.lineEdit_PSO_sigma->text();
	double sig = text3.toDouble();
	const QString text4  = this->ui.lineEdit_PSO_TolCon->text();
	double volPen = text4.toDouble();
	const QString text5  = this->ui.lineEdit_PSO_pop_size->text();
	int pop_size = text5.toInt();// number of particles (population size)
	const QString text6  = this->ui.lineEdit_PSO_max_velocity->text();
	double Vmax = text6.toDouble();// maximum velocity of a particle
	const QString text7  = this->ui.lineEdit_PSO_c1->text();
	double c1 = text7.toDouble();// cognition learning rate
	const QString text8  = this->ui.lineEdit_PSO_c2->text();
	double c2 = text8.toDouble();// social learning rate
	const QString text9  = this->ui.lineEdit_PSO_w->text();
	double w = text9.toDouble();// inertia weight
    //double volPen = 1e3;		// Penalty factor for volume constraint violation
	//////////

//////////////////////
			if(maxLfit<this->s->srf_m->b->L_max)
		{
			this->s->srf_m->flush_after_L(maxLfit);
		}
		else
		{
			if(maxLfit>this->s->srf_m->b->L_max)
			{
				this->set_L_max(maxLfit);		// in this case we need to actually increase the basis
			}
		}
		//this->s->srf_u->flush_after_L(maxLfit);
		if(this->gdim<this->gdim_fit_min)
		{
			int D_new = this->gdim_fit_min;
			//std::ostringstream os;
			os<<"Increasing Gaussian base point mesh to "<<D_new<<"  for accurate fitting. Please wait ..... ";report(os.str());
			this->s->set_new_gdim(D_new);
			this->gdim = D_new;
			os<<"Done!";report(os.str());
		}
//////////////////////////


	os<<"Maximum number of iterations: "<<niter<<"\nConstraint penalty factor: "<<volPen<<"\n";
	os<<"\nFitting coefficients up to Lmax: "<<maxLfit<<"\n"<<"\nPopulation size: "<<pop_size<<"\n"<<"------------------------------";report(os.str());os.str("");
	this->MCsig.resize(this->s->srf_m->xc.size());
	for(int i = 0;i<this->s->srf_m->xc.size();i++){this->MCsig[i] = sig;}	// fill MCsig with sig values
	
	// prepare the 
	double Eo = this->s->get_energy(1);	// this is the energy at the beginning (currently accepted)
	double Vo = this->s->srf_u->V;		// set the volume constraint to that of the starting configuration

	
	///// generate a starting set of particles (configurations) and initialize quantities
	std::vector<double> Energy_vec;	Energy_vec.reserve(niter);
	stop_optimization = false;
	double old_val = 0.0;		// here we will store the coefficient value before change
	int coord = 0;				// which coordinate was changed x,y,z
	int cix = 0;				// which coefficient was changed (linear index)
	int np = 0;					// particle index for particle generation
	int gbestix = 0;			// particle index of global best particle
	double gbestE = 1e20; 		// value of global best energy

	Eigen::MatrixXd  cwts;		// weights matrix for clks
	Eigen::MatrixXd xc_mx, yc_mx, zc_mx, E_mx;	// current configurations of particles and current energies
	Eigen::MatrixXd bxc_mx, byc_mx, bzc_mx;		// best configurations for all particles
	Eigen::MatrixXd v_xc, v_yc, v_zc;			// velocity vector for all particles for x y and z clks
	Eigen::MatrixXd xc_tmp, yc_tmp, zc_tmp;     // will hold the values of the starting shape

	xc_mx.resize(s->srf_m->xc.rows(),pop_size);xc_mx.fill(0.0);
	yc_mx.resize(s->srf_m->yc.rows(),pop_size);yc_mx.fill(0.0);
	zc_mx.resize(s->srf_m->zc.rows(),pop_size);zc_mx.fill(0.0);
	bxc_mx.resize(s->srf_m->xc.rows(),pop_size);bxc_mx.fill(0.0);
	byc_mx.resize(s->srf_m->yc.rows(),pop_size);byc_mx.fill(0.0);
	bzc_mx.resize(s->srf_m->zc.rows(),pop_size);bzc_mx.fill(0.0);
	v_xc.resize(s->srf_m->xc.rows(), pop_size);v_xc.fill(0.0);
	v_yc.resize(s->srf_m->yc.rows(), pop_size);v_yc.fill(0.0);
	v_zc.resize(s->srf_m->zc.rows(), pop_size);v_zc.fill(0.0);
	E_mx.resize(pop_size,1);E_mx.fill(1e40);

	xc_tmp = this->s->srf_m->xc;
	yc_tmp = this->s->srf_m->yc;
	zc_tmp = this->s->srf_m->zc;

	cwts.resize(s->srf_m->xc.rows(),1);cwts.fill(0.0);
	int ix = 0;
	for(int L = 0;L<=maxLfit;L++)						// loop over the L values up to maxLfit
	{
		for(int K = -L;K<=L;K++)						// loop over K
		{
            cwts(ix,0) = N_LK_bosh(L,K);
			ix++;
		}
	}

	///////////////////// PARTICLE INITIALIZATION LOOP //////////////////////
	bool particle_is_legal = false;
	while(np<pop_size & stop_optimization==0)
	{
		os<<"Initializing particle number: "<<np+1<<" ...";report(os.str());os.str("");
		particle_is_legal = false;
		while(!particle_is_legal & stop_optimization==0) //generate particles until we can accept one (i.e. does not intersect the vitteline membrane or itself)
		{
			QApplication::processEvents(QEventLoop::AllEvents);		//sosi -- in the future start a separate thread for the optimization
			//os<<"trying to initalize particle no.:\t"<<np+1<<"\t";report(os.str());os.str("");
			int maxix = (this->maxLfit + 1)*(this->maxLfit + 1);
			coord = int(unifRand(3));// generate a random integer between 1 and 3 to decide on the coordinate that will be changed x y or z
			cix = int(unifRand(maxix-1)); // random integer to choose the coefficient -- we are avoiding change of the first coefficient
			double new_val = 0.0;
			if(coord==1)
			{
				old_val = this->s->srf_m->xc(cix,0);
				new_val = randn_notrig(this->s->srf_m->xc(cix,0), this->MCsig[cix]);
				this->s->srf_m->xc(cix,0) = new_val;
			}
			if(coord==2)
			{
				old_val = this->s->srf_m->yc(cix,0);
				new_val = randn_notrig(this->s->srf_m->yc(cix,0), this->MCsig[cix]);
				this->s->srf_m->yc(cix,0) = new_val;
			}
			if(coord==3)
			{
				old_val = this->s->srf_m->zc(cix,0);
				new_val = randn_notrig(this->s->srf_m->zc(cix,0), this->MCsig[cix]);
				this->s->srf_m->zc(cix,0) = new_val;
			}
			this->s->update();
			 
			if(this->vtk_intersection_tests()==0)
			{
				os<<" Done!\n";report(os.str());os.str("");
				particle_is_legal = true;
				// store the configurations
				xc_mx.col(np) = this->s->srf_m->xc;
				yc_mx.col(np) = this->s->srf_m->yc;
				zc_mx.col(np) = this->s->srf_m->zc;
			}
			// restore working object to initial configuration
			if(coord==1){this->s->srf_m->xc(cix,0) = old_val;}
			if(coord==2){this->s->srf_m->yc(cix,0) = old_val;}
			if(coord==3){this->s->srf_m->zc(cix,0) = old_val;}
		}
		np++;	// increment the particle index counter
	}
	stop_optimization=false;
	///// the actual PSO iterations start here
	/////////////////////////////////////////
	double LARGE2 = 1e10;
	int iter = 0;
	bool volCon = 0;
	double Evol = 0.0;
	double E = 0.0;							// current energy
	this->update_vis_optim_hi_res();		// let's look at it before the iterations start

	
	////////////////////////////// START ITERATIONS //////////////////////////////////
	while(iter<niter & stop_optimization==0 )
	{
		os<<"Iteration\t"<<iter+1<<"\n";report(os.str());os.str("");
		///////////////////////// LOOP 1 --- calculate energies and determine best configurations
				np = 0;
				while(np<pop_size & stop_optimization==0)
				{
					QApplication::processEvents(QEventLoop::AllEvents);		//sosi -- in the future start a separate thread for the optimization
					// update the working object
					this->s->srf_m->set_clks(xc_mx.col(np), yc_mx.col(np), zc_mx.col(np));
					// calculate its energy
					E = this->s->get_energy() + volPen * std::abs(Vo-this->s->srf_m->V);
					os<<"Current energy: "<<E<<"\tCurrent best:  "<<gbestE;report(os.str());os.str("");
					this->update_vis_optim_hi_res();
					// check wether we need to make this configuration the best configuration for this particle
					if(E<E_mx(np,0) & this->vtk_intersection_tests() == 0)
					{
						E_mx(np,0) = E;
						bxc_mx.col(np) = xc_mx.col(np);
						byc_mx.col(np) = yc_mx.col(np);
						bzc_mx.col(np) = zc_mx.col(np);
						// check whether this configuration is the best in history
						if(E_mx(np,0)<gbestE)
						{
							gbestE = E_mx(np,0);
							gbestix =np;
							/// update the plots
							Energy_vec.push_back(gbestE);
							this->update_vis_optim_hi_res();//os<<" --- probability accept"<<fac1<<"\tProbability: "<<p;
							this->generate_XYplot(Energy_vec);
						}
					}
					np++;
				}
		///////////////////////// LOOP 2 --- calculate new particle positions
				np = 0;
				while(np<pop_size & stop_optimization==0)
				{
					QApplication::processEvents(QEventLoop::AllEvents);		//sosi -- in the future start a separate thread for the optimization
					// calculate particle velocities and new positions --- then the corresponding energy for every particle
					double rand1 = double(unifRand());
					double rand2 = double(unifRand());
					// update velocity
					v_xc.col(np) = w*v_xc.col(np).array() + c1 * Eigen::MatrixXd::Random(v_xc.rows(),1).array() * (bxc_mx.col(np)-xc_mx.col(np)).array() + c2 * Eigen::MatrixXd::Random(v_xc.rows(),1).array() * (bxc_mx.col(gbestix)-xc_mx.col(np)).array();
					v_yc.col(np) = w*v_yc.col(np).array() + c1 * Eigen::MatrixXd::Random(v_xc.rows(),1).array() * (byc_mx.col(np)-yc_mx.col(np)).array() + c2 * Eigen::MatrixXd::Random(v_xc.rows(),1).array() * (byc_mx.col(gbestix)-yc_mx.col(np)).array();
					v_zc.col(np) = w*v_zc.col(np).array() + c1 * Eigen::MatrixXd::Random(v_xc.rows(),1).array() * (bzc_mx.col(np)-zc_mx.col(np)).array() + c2 * Eigen::MatrixXd::Random(v_xc.rows(),1).array() * (bzc_mx.col(gbestix)-zc_mx.col(np)).array();
					// set upper limit for velocities
					v_xc = (v_xc.array()>Vmax).select(Vmax,v_xc.array());
					v_yc = (v_yc.array()>Vmax).select(Vmax,v_yc.array());
					v_zc = (v_zc.array()>Vmax).select(Vmax,v_zc.array());
					// update position
					xc_mx.col(np) = xc_mx.col(np) + v_xc.col(np);
					yc_mx.col(np) = yc_mx.col(np) + v_yc.col(np);
					zc_mx.col(np) = zc_mx.col(np) + v_zc.col(np);
					np++;
				}
		iter++;
	}
	os<<"Terminating particle swarm optimization -- updating to current best configuration\n";report(os.str());os.str("");
	this->s->srf_m->set_clks(xc_mx.col(gbestix), yc_mx.col(gbestix), zc_mx.col(gbestix));
	this->update_vis_optim_hi_res();//os<<" --- probability accept"<<fac1<<"\tProbability: "<<p;
	this->generate_XYplot(Energy_vec);
	
	/**/
}



/////////////////////// Objective and constraints functions for NLopt optimization routines ///////////////////////
double SimGast::myfunc(unsigned n, const double *x, double *grad, void *my_func_data)
{
	SimGast *app = (SimGast*) my_func_data;
	// check whether we need to stop the optimization or not
	QApplication::processEvents(QEventLoop::AllEvents);
	if(app->stop_optimization==1)
	{
		nlopt_force_stop(app->opt);
	}

	
	for(int i = 0;i<n;i++){app->X_o[app->parmix[i]] = x[i];} // set up the current full parameter vector
	int nc = app->X_o.size()/3;
	for(int i = 0;i<nc;i++)
	{
		app->s->srf_m->xc(i,0) = app->X_o[i] * app->tnlk[i];
		app->s->srf_m->yc(i,0) = app->X_o[i+nc] * app->tnlk[i];
		app->s->srf_m->zc(i,0) = app->X_o[i+2*nc] * app->tnlk[i];
	}
	app->s->srf_m->update_fast();
	double E = app->s->shell_energy();// calculate the energy
	//double E = 0.0;
	app->Energy_vec.push_back(E);
	// visualize
	if(app->func_counter % app->opt_plot_step == 0)		// don't plot every time the objective functionn is called
	{
		if(app->show_hi_res){app->update_vis_optim_hi_res();}
		else{app->update_vis_optim();}
		app->generate_XYplot(app->Energy_vec);
	}
	app->func_counter++;
	
	return E;

}
double SimGast::myfunc_unconstrained(unsigned n, const double *x, double *grad, void *my_func_data)
{
	// note: *x containts the fitting parameters only, this may not be the complete set of shape parameters
	//       therefore we reconstruct app->X_o from *x using app->parmix and then insert the current values
	//		 into app->s->srf_m. Also n contains the number of parameters in *x.

	SimGast *app = (SimGast*) my_func_data;
	// check whether we need to stop the optimization or not
	QApplication::processEvents(QEventLoop::AllEvents);
	if(app->stop_optimization==1)
	{
		nlopt_force_stop(app->opt);
	}

	
	for(int i = 0;i<n;i++){app->X_o[app->parmix[i]] = x[i];} // set up the current full parameter vector
	int nc = app->X_o.size()/3;
	for(int i = 0;i<nc;i++)
	{
		app->s->srf_m->xc(i,0) = app->X_o[i] * app->tnlk[i];
		app->s->srf_m->yc(i,0) = app->X_o[i+nc] * app->tnlk[i];
		app->s->srf_m->zc(i,0) = app->X_o[i+2*nc] * app->tnlk[i];
	}

	//app->s->srf_m->set_clks(xc, yc, zc);
	//app->s->srf_m->center_to_zero();
	app->s->update();
	double E = app->s->get_energy(0);
	app->Energy_vec.push_back(E);
	// add the contraints
	double f = 0.0;
	//double Ev = app->gamma_v * (1- (app->s->srf_m->V/app->s->srf_u->V))*(1- (app->s->srf_m->V/app->s->srf_u->V)); // original implementation
	double Ev = app->gamma_v * (app->s->srf_m->V- ((1+f)* app->s->srf_u->V))* (app->s->srf_m->V- ((1+f)* app->s->srf_u->V)); // Conte et al. 2008 (Equation 4)
	int intersection = (app->vtk_intersection_tests());
	double Einter = app->gamma_inter * double(intersection);		// energy in case of intersection
	double ET     = app->gamma_T * (1-app->s->srf_m->T) * (1-app->s->srf_m->T);// deviation of T from one is not allowed
	// visualize / report
	std::ostringstream os;
	//os<<"Vol. constraint: "<<std::abs((1-(app->s->srf_m->V/app->s->srf_u->V)))<<"\tIntersection : "<<intersection<<"\tT: "<<app->s->srf_m->T<<std::endl;app->report(os.str());os.str("");
	os<<"E:"<<E<<" -- Ev: "<<Ev<<"\tIntersec.: "<<intersection<<"\tT: "<<app->s->srf_m->T<<std::endl;app->report(os.str());os.str("");
	if(app->func_counter % app->opt_plot_step ==0)
	{
		if(app->show_hi_res){app->update_vis_optim_hi_res();}
		else{app->update_vis_optim();}
		app->generate_XYplot(app->Energy_vec);
	}
	app->func_counter++;
	// calculate the energy
	E +=(Ev+Einter + ET);
	return E;

}
typedef struct {
    int a;
	int nc;
	double Vo;
	SimGast *app;
} my_constraint_data;
double SimGast::myconstraint(unsigned n, const double *x, double *grad, void *data)
{	// Note gamma_T (check for consistency of shape calculation) has not been implemented here yet.
    my_constraint_data *d = (my_constraint_data *) data;
	SimGast *app = d->app;
	int a = d->a;
	//int nc = d->nc;
	double Vo = d->Vo;


	for(int i = 0;i<n;i++){app->X_o[app->parmix[i]] = x[i];} // set up the current full parameter vector
	int nc = app->X_o.size()/3;
	for(int i = 0;i<nc;i++)
	{
		app->s->srf_m->xc(i,0) = app->X_o[i] * app->tnlk[i];
		app->s->srf_m->yc(i,0) = app->X_o[i+nc] * app->tnlk[i];
		app->s->srf_m->zc(i,0) = app->X_o[i+2*nc] * app->tnlk[i];
	}
	//app->s->srf_m->set_clks(xc, yc, zc);
	//app->s->srf_m->center_to_zero();
	// calculate the constraint



	double cval =  0.0;	// constraint value
	std::ostringstream os;
	if(a==0)	// calculate volume and return volume constraint
	{
		app->s->srf_m->update_fast_volume();
		cval = 1- (app->s->srf_m->V/Vo);
		//os<<"Volume constraint: "<<cval<<"\tV diff.: "<<Vo-app->s->srf_m->V<<std::endl;app->report(os.str());os.str("");
	}
	if(a==1)	// check whether there is intersection (of any type) or not)
	{	
		
		app->s->update_fast();
		if(app->vtk_intersection_tests()>0) {cval = 1e20;}	// intersection test
		//os<<"Intersection: "<<cval<<std::endl;app->report(os.str());os.str("");
	}
    /*double a = d->a, b = d->b;
    if (grad) {
        grad[0] = 3 * a * (a*x[0] + b) * (a*x[0] + b);
        grad[1] = -1.0;
    }
    return ((a*x[0] + b) * (a*x[0] + b) * (a*x[0] + b) - x[1]);
	*/
	return cval;
 }
//////////////////////// FITTING --- COBYLA ////////////////
void SimGast::on_pushButton_start_COBYLA_clicked()
{	
		// prepare the shell morphology for the limited Lmax specified for fitting
	
		std::ostringstream os;os<<"--------- STARTING COBYLA OPTIMIZATION --------------"<<std::endl;s->disp(os,*s);;//report(os.str());
		stop_optimization=false;
		// prepare the shell morphology for the Lmax specified for fitting and make sure it has minimum gdim
		const QString text1  = this->ui.lineEdit_COBYLA_Lmax_fit->text();
		maxLfit = text1.toInt();
		if(maxLfit<this->s->srf_m->b->L_max)
		{
			this->s->srf_m->flush_after_L(maxLfit);
		}
		else
		{
			if(maxLfit>this->s->srf_m->b->L_max)
			{
				this->set_L_max(maxLfit);		// in this case we need to actually increase the basis
			}
		}
		//this->s->srf_u->flush_after_L(maxLfit);
		if(this->gdim<this->gdim_fit_min)
		{
			int D_new = this->gdim_fit_min;
			std::ostringstream os;os<<"Increasing Gaussian base point mesh to "<<D_new<<"  for accurate fitting. Please wait ..... ";report(os.str());
			this->s->set_new_gdim(D_new);
			this->gdim = D_new;
			os<<"Done!";report(os.str());
		}



		int old_tri_n = this->s->srf_m->sm->tri_n;
		this->s->srf_m->set_new_spherical_mesh(tri_n_intersection_test);		// we use a separate resolution mesh for intersection tests
		this->s->srf_m->center_to_zero();
		double E = this->s->get_energy(1);
		//this->s->prepare_energy_calc();
		//this->s->update();
		os<<"\n"<<"Active region area fraction:"<<s->active_area_fraction<<std::endl;
		//os<<"Vout Vin Aout Ain:"<<<<"\n"<<s->vin<<"\t"<<s->aout<<"\t"<<s->ain<<"\t"<<std::endl;
		os<<"Volume enclosed by outer surface\t"<<s->vout<<std::endl;
		os<<"Area of outer surface\t"<<s->aout<<std::endl;
		os<<"Volume enclosed by inner surface\t\t"<<s->vin<<std::endl;
		os<<"Area of inner surface\t"<<s->ain<<std::endl;

		os<<"Shear and stretch: "<<this->s->E_nHook<<"\tBending: "<<this->s->Eb<<"\nTotal energy:\n"<<E<<std::endl;
		report(os.str());os.str("");

		this->COBYLA_start();
		this->s->update();
		this->on_pushButton_status_clicked();
		if(stop_optimization){os<<"Optimization stopped at user's request";report(os.str());os.str("");}
		this->s->srf_m->set_new_spherical_mesh(old_tri_n);		// reset the spherical mesh to what it was before
}


void SimGast::COBYLA_start()
{	
	std::ostringstream os;//os<<"---------STARTING COBYLA OPTIMIZATION-----------"<<std::endl;
	Energy_vec.clear();
	func_counter = 0;
	// define preliminary quantities
	const QString text1  = this->ui.lineEdit_COBYLA_max_fun->text();
	int maxfun = text1.toInt();
	const QString text2  = this->ui.lineEdit_COBYLA_xtol->text();
	double xtol = text2.toDouble();
	const QString text3  = this->ui.lineEdit_COBYLA_dx->text();
	double dxo = text3.toDouble();
	const QString text4  = this->ui.lineEdit_COBYLA_vtol->text();
	double vtol = text4.toDouble();
	double Eo = this->s->get_energy(0);	// this is the energy at the beginning (currently accepted)
	double Vo = this->s->srf_u->V;		// set the volume constraint to that of the starting configuration
	int nparms = 3*(maxLfit + 1)*(maxLfit + 1);
	int nc = (maxLfit + 1)*(maxLfit + 1);
	Energy_vec.push_back(Eo);
	// define function pointers
	double (*fptr) (unsigned n, const double *x, double *grad, void *my_func_data) = &SimGast::myfunc;
	double (*cptr) (unsigned n, const double *x, double *grad, void *data)			= &SimGast::myconstraint;

		// specify a starting vector
	onlk.resize(nc);over_NLK(maxLfit, onlk);
	tnlk.resize(nc);times_NLK(maxLfit, tnlk);

	// sosi
	std::vector<double> tmp1(nc,1.0);
	std::vector<double> tmp2(nc,1.0);
	onlk = tmp1;
	tnlk = tmp2;

	// generate the fitting parameter index matrices
	this->set_L_max(maxLfit);
	if(this->ui.checkBox_symmetry->isChecked()){this->enforce_symmetry();};
	X_o.resize(0);
	X_o.reserve(nparms);//std::fill(X_o.begin(), X_o.end(), 0);
	for(int i = 0;i<nc;i++){X_o.push_back(this->s->srf_m->xc(i,0)*onlk[i]);}
	for(int i = 0;i<nc;i++){X_o.push_back(this->s->srf_m->yc(i,0)*onlk[i]);}
	for(int i = 0;i<nc;i++){X_o.push_back(this->s->srf_m->zc(i,0)*onlk[i]);}
	parmix.resize(0);
	parmix.reserve(nparms);
	///
	int counter = 0;
	std::vector<double> stvec;
	stvec.reserve(nparms);
	for(int i = 0;i<nc;i++)
	{
		if(ixc[i])	
		{
			stvec.push_back(this->s->srf_m->xc(i,0)*onlk[i]);
			parmix.push_back(counter);// store the index (into X_o) of the added parameter
			
		}	
		counter++;
	}						
	
	for(int i = 0;i<nc;i++){if(iyc[i])	{stvec.push_back(this->s->srf_m->yc(i,0)*onlk[i]);parmix.push_back(counter);}counter++;	}
	for(int i = 0;i<nc;i++){if(izc[i])	{stvec.push_back(this->s->srf_m->zc(i,0)*onlk[i]);parmix.push_back(counter);}counter++;	}
	os<<"Number of fitting parameters: "<< stvec.size()<< " of possible "<< this->s->srf_m->xc.rows()*3<<"\n";report(os.str());os.str("");




	// configure the optimization object
	//nlopt_opt opt;
	opt = nlopt_create(NLOPT_LN_COBYLA, stvec.size()); // algorithm and dimensionality 
	//opt = nlopt_create(NLOPT_LN_BOBYQA, nparms); // algorithm and dimensionality 
	//opt = nlopt_create(NLOPT_LN_NELDERMEAD, nparms); // algorithm and dimensionality
	nlopt_set_min_objective(opt, *fptr, this);		//specify the objective function// pass the "this" pointer as the data parameter (gets cast in myfunc)
	my_constraint_data data[2] = { {0,nc,Vo, this}, {1,nc,Vo, this} };
	nlopt_add_equality_constraint(opt, *cptr, &data[0], vtol);
	nlopt_add_equality_constraint(opt, *cptr, &data[1], vtol);

	// specify additional parameters
	std::vector<double> dx(stvec.size(),dxo);
	nlopt_set_initial_step(opt, &dx[0]);
	nlopt_set_xtol_rel(opt, xtol);
	nlopt_set_maxeval(opt, maxfun);


	double minf; // the minimum objective value, upon return

	// test the objective function
	minf = fptr(stvec.size(), &stvec[0], NULL, this);

	// start the optimization
	clock_t begin,end;
	double time_spent;
	begin = clock();

	if (nlopt_optimize(opt, &stvec[0], &minf) < 0) {
		os<<"COBYLA stopped\n";report(os.str());os.str("");
	}
	else {
		os<<"COBYLA found minimum!"<<minf<<std::endl;report(os.str());os.str("");
		//printf("found minimum at f(%g,%g) = %0.10g\n", x[0], x[1], minf);
	}
	end = clock();
	time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
	os<<"Time spent: "<<time_spent<<" seconds"<<std::endl;report(os.str());os.str("");

	minf = fptr(stvec.size(), &stvec[0], NULL, this);		// call one last time with the final result
	nlopt_destroy(opt);
	stvec.clear();
	dx.clear();
	/**/
}


//////////////////////// FITTING --- Sbplx
void SimGast::on_pushButton_start_SBPLX_clicked()
{	
		std::ostringstream os;os<<"--------- STARTING SUBPLEX FITTING PROCEDURE --------------"<<std::endl;s->disp(os,*s);//report(os.str());
		// check whether we have fields defined or not
		bool fields_defined =  (s->sfindx[0] >-1 | s->sfindx[1] >-1 | s->sfindx[2] >-1 | s->sfindx[3] >-1 | s->sfindx[4] >-1 | s->sfindx[5] >-1); // i.e. if any of them points to an existing field

		if(fields_defined)
		{
		stop_optimization=false;
		// prepare the shell morphology for the Lmax specified for fitting and make sure it has minimum gdim
		const QString text1  = this->ui.lineEdit_SBPLX_Lmax_fit->text();
		maxLfit = text1.toInt();
		if(maxLfit<this->s->srf_m->b->L_max)
		{
			this->s->srf_m->flush_after_L(maxLfit);
		}
		else
		{
			if(maxLfit>this->s->srf_m->b->L_max)
			{
				this->set_L_max(maxLfit);		// in this case we need to actually increase the basis
			}
		}
		//this->s->srf_u->flush_after_L(maxLfit);
		if(this->gdim<this->gdim_fit_min)
		{
			int D_new = this->gdim_fit_min;
			std::ostringstream os;os<<"Increasing Gaussian base point mesh to "<<D_new<<"  for accurate fitting. Please wait ..... ";report(os.str());
			this->s->set_new_gdim(D_new);
			this->gdim = D_new;
			os<<"Done!";report(os.str());
		}

		int old_tri_n = this->s->srf_m->sm->tri_n;
		this->s->srf_m->set_new_spherical_mesh(tri_n_intersection_test);		// we use a low resolution mesh for intersection tests
		this->s->srf_m->center_to_zero();
		double E = this->s->get_energy(1);
		os<<"\n"<<"Active region area fraction:"<<s->active_area_fraction<<std::endl;
		//os<<"Vout Vin Aout Ain:"<<<<"\n"<<s->vin<<"\t"<<s->aout<<"\t"<<s->ain<<"\t"<<std::endl;
		os<<"Volume enclosed by outer surface\t"<<s->vout<<std::endl;
		os<<"Area of outer surface\t"<<s->aout<<std::endl;
		os<<"Volume enclosed by inner surface\t\t"<<s->vin<<std::endl;
		os<<"Area of inner surface\t"<<s->ain<<std::endl;

		os<<"Shear and stretch: "<<this->s->E_nHook<<"\tBending: "<<this->s->Eb<<"\nTotal energy:\n"<<E<<std::endl;
		report(os.str());os.str("");

		this->SBPLX_start();
		this->s->update_scalar_fields();
		this->s->update();
		this->on_pushButton_status_clicked();
		if(stop_optimization){os<<"Optimization stopped at user's request";report(os.str());os.str("");}
		this->s->srf_m->set_new_spherical_mesh(old_tri_n);		// reset the spherical mesh to what it was before
		this->update_vtk();

		}
}


void SimGast::SBPLX_start()
{	
	std::ostringstream os;//os<<"---------STARTING SBPLX OPTIMIZATION-----------"<<std::endl;
	Energy_vec.clear();
	func_counter = 0;
	// define preliminary quantities
	const QString text1  = this->ui.lineEdit_SBPLX_max_fun->text();
	int maxfun = text1.toInt();
	const QString text2  = this->ui.lineEdit_SBPLX_gamma_inter->text();
	gamma_inter = text2.toDouble();
	const QString text3  = this->ui.lineEdit_SBPLX_dx->text();
	double dxo = text3.toDouble();
	const QString text4  = this->ui.lineEdit_SBPLX_gamma_v->text();
	gamma_v = text4.toDouble();
	double Eo = this->s->get_energy(0);	// this is the energy at the beginning (currently accepted)
	double Vo = this->s->srf_u->V;		// set the volume constraint to that of the starting configuration
	int nparms = 3*(maxLfit + 1)*(maxLfit + 1);
	int nc = (maxLfit + 1)*(maxLfit + 1);
	Energy_vec.push_back(Eo);
	// define function pointers
	double (*fptr) (unsigned n, const double *x, double *grad, void *my_func_data) = &SimGast::myfunc_unconstrained;
	//double (*cptr) (unsigned n, const double *x, double *grad, void *data)			= &SimGast::myconstraint;
	
	// specify a starting vector
	onlk.resize(nc);over_NLK(maxLfit, onlk);
	tnlk.resize(nc);times_NLK(maxLfit, tnlk);
	//sosi
	std::vector<double> tmp1(nc,1.0);
	std::vector<double> tmp2(nc,1.0);
	onlk = tmp1;
	tnlk = tmp2;
	
	// generate fitting parameter index matrices
	this->set_L_max(maxLfit);
	//if(this->ui.checkBox_symmetry->isChecked()){this->enforce_symmetry();};
	X_o.resize(0);
	X_o.reserve(nparms);//std::fill(X_o.begin(), X_o.end(), 0);
	for(int i = 0;i<nc;i++){X_o.push_back(this->s->srf_m->xc(i,0)*onlk[i]);}
	for(int i = 0;i<nc;i++){X_o.push_back(this->s->srf_m->yc(i,0)*onlk[i]);}
	for(int i = 0;i<nc;i++){X_o.push_back(this->s->srf_m->zc(i,0)*onlk[i]);}
	parmix.resize(0);
	parmix.reserve(nparms);
	///
	int counter = 0;
	std::vector<double> stvec;
	stvec.reserve(nparms);
	for(int i = 0;i<nc;i++)
	{
		if(ixc[i])	
		{
			stvec.push_back(this->s->srf_m->xc(i,0)*onlk[i]);
			parmix.push_back(counter);// store the index (into X_o) of the added parameter
			
		}	
		counter++;
	}						
	
	for(int i = 0;i<nc;i++){if(iyc[i])	{stvec.push_back(this->s->srf_m->yc(i,0)*onlk[i]);parmix.push_back(counter);}counter++;	}
	for(int i = 0;i<nc;i++){if(izc[i])	{stvec.push_back(this->s->srf_m->zc(i,0)*onlk[i]);parmix.push_back(counter);}counter++;	}
	os<<"Number of fitting parameters: "<< stvec.size()<< " of possible "<< this->s->srf_m->xc.rows()*3<<"\n";report(os.str());os.str("");

	// configure the optimization object
	opt = nlopt_create(NLOPT_LN_SBPLX, stvec.size()); // algorithm and dimensionality 
	//opt = nlopt_create(NLOPT_LN_BOBYQA, nparms); // algorithm and dimensionality 
	//opt = nlopt_create(NLOPT_LN_NELDERMEAD, nparms); // algorithm and dimensionality
	nlopt_set_min_objective(opt, *fptr, this);		//specify the objective function// pass the "this" pointer as the data parameter (gets cast in myfunc)
	//my_constraint_data data[2] = { {0,nc,Vo, this}, {1,nc,Vo, this} };
	//nlopt_add_equality_constraint(opt, *cptr, &data[0], vtol);
	//nlopt_add_equality_constraint(opt, *cptr, &data[1], vtol);

	// specify additional parameters
	std::vector<double> dx(stvec.size(),dxo);
	nlopt_set_initial_step(opt, &dx[0]);
	//nlopt_set_xtol_rel(opt, xtol);
	nlopt_set_maxeval(opt, maxfun);



	double minf; // the minimum objective value, upon return

	// test the objective function
	//minf = fptr(nparms, &stvec[0], NULL, this);

	// start the optimization
	clock_t begin,end;
	double time_spent;
	begin = clock();

	if (nlopt_optimize(opt, &stvec[0], &minf) < 0) {
		os<<"SBPLX stopped\n";report(os.str());os.str("");
	}
	else {
		os<<"SBPLX found minimum!"<<minf<<std::endl;report(os.str());os.str("");
		//printf("found minimum at f(%g,%g) = %0.10g\n", x[0], x[1], minf);
	}
	end = clock();
	time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
	os<<"Time spent: "<<time_spent<<" seconds"<<std::endl;report(os.str());os.str("");

	minf = fptr(stvec.size(), &stvec[0], NULL, this);		// call one last time with the final result

	nlopt_destroy(opt);
	stvec.clear();
	dx.clear();
	/**/
}


//////////////////////// FITTING --- ISRES
void SimGast::on_pushButton_start_ISRES_clicked()
{	
	stop_optimization=false;
	if(this->ui.comboBox_sf1->currentIndex()>0 & this->ui.comboBox_sf2->currentIndex()>0 & this->ui.comboBox_sf3->currentIndex()>0)
	{
		// prepare the shell morphology for the limited Lmax specified for fitting
		const QString text1  = this->ui.lineEdit_ISRES_Lmax_fit->text();
		maxLfit = text1.toInt();
		this->s->srf_m->flush_after_L(maxLfit);
		//this->s->srf_u->flush_after_L(maxLfit);
		
		std::ostringstream os;os<<"---------STARTING CONFIGURATION --------------"<<std::endl;s->disp(os,*s);//report(os.str());
		int old_tri_n = this->s->srf_m->sm->tri_n;
		this->s->srf_m->set_new_spherical_mesh(tri_n_intersection_test);		// we use a low resolution mesh for intersection tests
		this->s->srf_m->center_to_zero();
		double E = this->s->get_energy(1);
		os<<"\n"<<"Active region area fraction:"<<s->active_area_fraction<<std::endl;
		//os<<"Vout Vin Aout Ain:"<<<<"\n"<<s->vin<<"\t"<<s->aout<<"\t"<<s->ain<<"\t"<<std::endl;
		os<<"Volume enclosed by outer surface\t"<<s->vout<<std::endl;
		os<<"Area of outer surface\t"<<s->aout<<std::endl;
		os<<"Volume enclosed by inner surface\t\t"<<s->vin<<std::endl;
		os<<"Area of inner surface\t"<<s->ain<<std::endl;

		os<<"Shear and stretch: "<<this->s->E_nHook<<"\tBending: "<<this->s->Eb<<"\nTotal energy:\n"<<E<<std::endl;
		report(os.str());os.str("");

		this->ISRES_start();
		this->s->update();
		this->on_pushButton_status_clicked();
		if(stop_optimization){os<<"Optimization stopped at user's request";report(os.str());os.str("");}
		this->s->srf_m->set_new_spherical_mesh(old_tri_n);		// reset the spherical mesh to what it was before
		this->update_vtk();
	}
	else
	{
		this->message("Please select three scalar fields first.");
	}

}


void SimGast::ISRES_start()
{	
	std::ostringstream os;os<<"---------STARTING ISRES OPTIMIZATION-----------"<<std::endl;
	Energy_vec.clear();
	func_counter = 0;
	// define preliminary quantities
	const QString text1  = this->ui.lineEdit_ISRES_max_fun->text();
	int maxfun = text1.toInt();
	const QString text2  = this->ui.lineEdit_ISRES_pop->text();
	int pop = text2.toInt();
	const QString text4  = this->ui.lineEdit_ISRES_vtol->text();
	double vtol = text4.toDouble();

	const QString text5  = this->ui.lineEdit_ISRES_lb->text();
	double lb_const = text5.toDouble();
	const QString text6  = this->ui.lineEdit_ISRES_ub->text();
	double ub_const = text6.toDouble();

	double Eo = this->s->get_energy(0);	// this is the energy at the beginning (currently accepted)
	double Vo = this->s->srf_u->V;		// set the volume constraint to that of the starting configuration
	int nparms = 3*(maxLfit + 1)*(maxLfit + 1);
	int nc = (maxLfit + 1)*(maxLfit + 1);
	Energy_vec.push_back(Eo);
	// define function pointers
	double (*fptr) (unsigned n, const double *x, double *grad, void *my_func_data) = &SimGast::myfunc;
	double (*cptr) (unsigned n, const double *x, double *grad, void *data)			= &SimGast::myconstraint;

	// specify a starting vector, lower and upper bounds
	onlk.resize(nc);over_NLK(maxLfit, onlk);
	tnlk.resize(nc);times_NLK(maxLfit, tnlk);
	// sosi
	std::vector<double> tmp1(nc,1.0);
	std::vector<double> tmp2(nc,1.0);
	onlk = tmp1;
	tnlk = tmp2;
	///
// generate the fitting parameter index matrices
	this->set_L_max(maxLfit);
	if(this->ui.checkBox_symmetry->isChecked()){this->enforce_symmetry();};
	X_o.resize(0);
	X_o.reserve(nparms);//std::fill(X_o.begin(), X_o.end(), 0);
	for(int i = 0;i<nc;i++){X_o.push_back(this->s->srf_m->xc(i,0)*onlk[i]);}
	for(int i = 0;i<nc;i++){X_o.push_back(this->s->srf_m->yc(i,0)*onlk[i]);}
	for(int i = 0;i<nc;i++){X_o.push_back(this->s->srf_m->zc(i,0)*onlk[i]);}
	parmix.resize(0);
	parmix.reserve(nparms);
	///
	int counter = 0;
	std::vector<double> stvec;
	stvec.reserve(nparms);
	std::vector<double> lb;lb.reserve(nparms);
	std::vector<double> ub;ub.reserve(nparms);
	for(int i = 0;i<nc;i++)
	{
		if(ixc[i])	
		{
			stvec.push_back(this->s->srf_m->xc(i,0)*onlk[i]);
			lb.push_back(this->s->srf_m->xc(i,0)*onlk[i] + lb_const);ub.push_back(this->s->srf_m->xc(i,0)*onlk[i] + ub_const);
			parmix.push_back(counter);// store the index (into X_o) of the added parameter
		}	
		counter++;
	}						
	
	for(int i = 0;i<nc;i++){if(iyc[i])	{stvec.push_back(this->s->srf_m->yc(i,0)*onlk[i]);parmix.push_back(counter);lb.push_back(this->s->srf_m->yc(i,0)*onlk[i] + lb_const);ub.push_back(this->s->srf_m->yc(i,0)*onlk[i] + ub_const);}counter++;	}
	for(int i = 0;i<nc;i++){if(izc[i])	{stvec.push_back(this->s->srf_m->zc(i,0)*onlk[i]);parmix.push_back(counter);lb.push_back(this->s->srf_m->zc(i,0)*onlk[i] + lb_const);ub.push_back(this->s->srf_m->zc(i,0)*onlk[i] + ub_const);}counter++;	}
	os<<"Number of fitting parameters: "<< stvec.size()<< " of possible "<< this->s->srf_m->xc.rows()*3<<"\n";report(os.str());os.str("");

	// configure the optimization object
	opt = nlopt_create(NLOPT_GN_ISRES, stvec.size()); // algorithm and dimensionality 
	//opt = nlopt_create(NLOPT_LN_BOBYQA, nparms); // algorithm and dimensionality 
	//opt = nlopt_create(NLOPT_LN_NELDERMEAD, nparms); // algorithm and dimensionality
	nlopt_set_min_objective(opt, *fptr, this);		//specify the objective function// pass the "this" pointer as the data parameter (gets cast in myfunc)
	my_constraint_data data[2] = { {0,nc,Vo, this}, {1,nc,Vo, this} };
	nlopt_add_equality_constraint(opt, *cptr, &data[0], vtol);
	nlopt_add_equality_constraint(opt, *cptr, &data[1], vtol);
	// specify additional parameters
	nlopt_set_lower_bounds(opt, &lb[0]);
	nlopt_set_upper_bounds(opt, &ub[0]);
	nlopt_set_population(opt, pop);
	nlopt_set_maxeval(opt, maxfun);

	double minf; // the minimum objective value, upon return

	// test the objective function
	minf = fptr(stvec.size(), &stvec[0], NULL, this);

	// start the optimization
	clock_t begin,end;
	double time_spent;
	begin = clock();

	if (nlopt_optimize(opt, &stvec[0], &minf) < 0) {
		os<<"ISRES stopped\n";report(os.str());os.str("");
	}
	else {
		os<<"ISRES found minimum!"<<minf<<std::endl;report(os.str());os.str("");
		//printf("found minimum at f(%g,%g) = %0.10g\n", x[0], x[1], minf);
	}
	end = clock();
	time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
	os<<"Time spent: "<<time_spent<<" seconds"<<std::endl;report(os.str());os.str("");


	minf = fptr(stvec.size(), &stvec[0], NULL, this);		// call one last time with the final result

	// cleanup
	nlopt_destroy(opt);
	lb.clear();
	ub.clear();
	stvec.clear();
	
	/**/
}

