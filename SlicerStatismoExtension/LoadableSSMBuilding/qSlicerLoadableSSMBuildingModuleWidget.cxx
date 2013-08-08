/*==============================================================================

  Program: 3D Slicer

  Portions (c) Copyright Brigham and Women's Hospital (BWH) All Rights Reserved.

  See COPYRIGHT.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.

==============================================================================*/

// Qt includes
#include <QDebug>
#include <QFileDialog>

// SlicerQt includes
#include "qSlicerLoadableSSMBuildingModuleWidget.h"
#include "ui_qSlicerLoadableSSMBuildingModuleWidget.h"

// Statismo includes
#include "statismo/StatisticalModel.h"
#include "Representers/VTK/vtkPolyDataRepresenter.h"
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>

#include "itkMeshRepresenter.h"
#include "itkMesh.h"
#include "statismo_ITK/itkStatisticalModel.h"
#include "itkMeshFileWriter.h"

//Library for GPFittingMesh
#include "itkCovarianceFunctionModelBuilder.h"
#include "statismo_ITK/itkDataManager.h"
#include <sys/types.h>
#include <errno.h>

//library for buildPartiallyFixedModel
#include "itkDirectory.h"
#include "itkPointsLocator.h"
#include "itkMeshFileReader.h"
#include "statismo_ITK/itkPartiallyFixedModelBuilder.h"

//library for modelFitting
#include "itkLinearInterpolateImageFunction.h"
#include "statismo_ITK/itkStatisticalShapeModelTransform.h"
#include "itkPenalizingMeanSquaresPointSetToImageMetric.h"
#include "itkLBFGSOptimizer.h"
#include "itkPointSetToPointSetRegistrationMethod.h"
#include "itkPointSetToImageRegistrationMethod.h"
#include "itkImageFileReader.h"
#include "itkCommand.h"


// LoadableSSMBuilding Logic includes
#include "vtkSlicerLoadableSSMBuildingLogic.h"

#include "vtkMRMLModelNode.h"
#include "vtkMRMLModelDisplayNode.h"
#include "vtkMRMLScene.h"
#include "vtkMRMLLayoutNode.h"
#include "vtkMRMLDoubleArrayNode.h"
#include <vtkDoubleArray.h>
#include "vtkMRMLChartViewNode.h"
#include "vtkMRMLChartNode.h"


using namespace statismo;
typedef vtkPolyDataRepresenter VtkRepresenterType;
typedef StatisticalModel<VtkRepresenterType> VtkStatisticalModelType;

const unsigned Dimensions = 3;
typedef itk::Mesh<float, Dimensions> MeshType;
typedef itk::MeshFileReader<MeshType> MeshFileReaderType;
typedef itk::MeshRepresenter<float, Dimensions> ItkRepresenterType;
typedef itk::StatisticalModel<ItkRepresenterType> ItkStatisticalModelType;
typedef vnl_vector<statismo::ScalarType> itkVectorType;


//-----------------------------------------------------------------------------
/// \ingroup Slicer_QtModules_ExtensionTemplate
class qSlicerLoadableSSMBuildingModuleWidgetPrivate: public Ui_qSlicerLoadableSSMBuildingModuleWidget
{
public:
  qSlicerLoadableSSMBuildingModuleWidgetPrivate();

  VtkStatisticalModelType* vtkStatModel;
  void setVtkStatModel(VtkStatisticalModelType* newVtkStatModel);
  VtkStatisticalModelType* getVtkStatModel();
  ItkStatisticalModelType* itkStatModel;
  void setItkStatModel(ItkStatisticalModelType* newItkStatModel);
  ItkStatisticalModelType* getItkStatModel();
};

//-----------------------------------------------------------------------------
// qSlicerLoadableSSMBuildingModuleWidgetPrivate methods

//-----------------------------------------------------------------------------
qSlicerLoadableSSMBuildingModuleWidgetPrivate::qSlicerLoadableSSMBuildingModuleWidgetPrivate()
{
}

void qSlicerLoadableSSMBuildingModuleWidgetPrivate::setVtkStatModel(VtkStatisticalModelType* newVtkStatModel)
{
  this->vtkStatModel = newVtkStatModel;
}

VtkStatisticalModelType* qSlicerLoadableSSMBuildingModuleWidgetPrivate::getVtkStatModel()
{
  VtkStatisticalModelType* newVtkStatModel;
  return newVtkStatModel = this->vtkStatModel;
}

void qSlicerLoadableSSMBuildingModuleWidgetPrivate::setItkStatModel(ItkStatisticalModelType* newItkStatModel)
{
  this->itkStatModel = newItkStatModel;
}

ItkStatisticalModelType* qSlicerLoadableSSMBuildingModuleWidgetPrivate::getItkStatModel()
{
  ItkStatisticalModelType* newItkStatModel;
  return newItkStatModel = this->itkStatModel;
}
//-----------------------------------------------------------------------------
// qSlicerLoadableSSMBuildingModuleWidget methods

//-----------------------------------------------------------------------------
qSlicerLoadableSSMBuildingModuleWidget::qSlicerLoadableSSMBuildingModuleWidget(QWidget* _parent)
  : Superclass( _parent )
  , d_ptr( new qSlicerLoadableSSMBuildingModuleWidgetPrivate )
{
}

//-----------------------------------------------------------------------------
qSlicerLoadableSSMBuildingModuleWidget::~qSlicerLoadableSSMBuildingModuleWidget()
{
}

//-----------------------------------------------------------------------------
void qSlicerLoadableSSMBuildingModuleWidget::setup()
{
  Q_D(qSlicerLoadableSSMBuildingModuleWidget);
  d->setupUi(this);
  
}

//-----------------------------------------------------------------------------
void qSlicerLoadableSSMBuildingModuleWidget::onSelectInputModel()
{
  Q_D(qSlicerLoadableSSMBuildingModuleWidget);
  QString inputFile = QFileDialog::getOpenFileName(this, "Select input model", QString());
  d->modelNamePath->setText(inputFile);
  // Load the model
  std::string modelString = inputFile.toStdString();
  if (d->radioButtonVTK->isChecked()){
    try {
      d->setVtkStatModel(d->vtkStatModel->Load(modelString.c_str()));

      // Display Eigen spectrum
      displayEigenSpectrum();
      // Set the number of components for the pcSlider
      d->pcSlider->setMaximum(d->vtkStatModel->GetNumberOfPrincipalComponents());

      //Calculate mean
      vtkPolyData* meanModel = d->vtkStatModel->DrawMean();
      // Add mean model to the scene
      vtkSlicerLoadableSSMBuildingLogic* moduleLogic = vtkSlicerLoadableSSMBuildingLogic::New();
      moduleLogic->DisplaySampleModel(meanModel, this->mrmlScene());

    }
	  catch (StatisticalModelException& e) {
		  std::cout << "Exception occured while building the shape model" << std::endl;
		  std::cout << e.what() << std::endl;
	  }
  }
  if (d->radioButtonITK->isChecked()){
    try {
      ItkStatisticalModelType* modelITK = ItkStatisticalModelType::New();
      modelITK->Load(modelString.c_str());
      d->setItkStatModel(modelITK);
      //d->setItkStatModel(d->itkStatModel->Load(modelString.c_str()));

      // Display Eigen spectrum
      displayEigenSpectrum();
      // Set the number of components for the pcSlider
      d->pcSlider->setMaximum(d->itkStatModel->GetNumberOfPrincipalComponents());

      //Calculate mean
      MeshType::Pointer meanItkModel = d->itkStatModel->DrawMean();
	     itk::MeshFileWriter<MeshType>::Pointer writer = itk::MeshFileWriter<MeshType>::New();
	     writer->SetFileName("meanItkModel.vtk");
	     writer->SetInput(meanItkModel);
	     writer->Update();

	     vtkPolyDataReader* reader = vtkPolyDataReader::New();
      reader->SetFileName("meanItkModel.vtk");
      reader->Update();
      vtkPolyData* meanModel = vtkPolyData::New();
      meanModel->ShallowCopy(reader->GetOutput());

      // Add mean model to the scene
      vtkSlicerLoadableSSMBuildingLogic* moduleLogic = vtkSlicerLoadableSSMBuildingLogic::New();
      moduleLogic->DisplaySampleModel(meanModel, this->mrmlScene());

    }
	  catch (itk::ExceptionObject& o) {
		  std::cout << "Exception occured while building the shape model" << std::endl;
		  std::cout << o << std::endl;
	 }
 }
  
}

//-----------------------------------------------------------------------------
void qSlicerLoadableSSMBuildingModuleWidget::onSelect()
{
  /**** VTK Model ****/
  using std::auto_ptr;
  Q_D(qSlicerLoadableSSMBuildingModuleWidget);
  vtkPolyData* samplePC = vtkPolyData::New();
  // Get the model name selected by the user
  if (d->radioButtonVTK->isChecked()){
    d->vtkStatModel = d->getVtkStatModel();
    int nbPrincipalComponent = d->vtkStatModel->GetNumberOfPrincipalComponents();
    VectorType coefficients = VectorType::Zero(nbPrincipalComponent);
    int pc = static_cast<int>(d->pcSlider->value())-1; // -1 because user can choose between 1 and max
    coefficients(pc) = d->stdSlider->value();
    samplePC = d->vtkStatModel->DrawSample(coefficients);
    double prob = d->vtkStatModel->ComputeProbabilityOfDataset(samplePC);
    std::cout<<"prob = "<<prob<<std::endl;
  }
  if (d->radioButtonITK->isChecked()){
    d->itkStatModel = d->getItkStatModel();
    int nbPrincipalComponent = d->itkStatModel->GetNumberOfPrincipalComponents();
    itkVectorType coefficients(nbPrincipalComponent,0.0); // set the vector to 0
    int pc = static_cast<int>(d->pcSlider->value())-1; // -1 because user can choose between 1 and max
    coefficients(pc) = d->stdSlider->value();
    MeshType::Pointer itkSamplePC = d->itkStatModel->DrawSample(coefficients);
    
    double prob = d->itkStatModel->ComputeProbabilityOfDataset(itkSamplePC);
    std::cout<<"prob = "<<prob<<std::endl;

    itk::MeshFileWriter<MeshType>::Pointer writer = itk::MeshFileWriter<MeshType>::New();
    writer->SetFileName("itkSamplePC.vtk");
    writer->SetInput(itkSamplePC);
    writer->Update();

    vtkPolyDataReader* reader = vtkPolyDataReader::New();
    reader->SetFileName("itkSamplePC.vtk");
    reader->Update();
    samplePC->ShallowCopy(reader->GetOutput());
  } 
  // Add polydata sample to the scene
  vtkSlicerLoadableSSMBuildingLogic* moduleLogic = vtkSlicerLoadableSSMBuildingLogic::New();
  moduleLogic->DisplaySampleModel(samplePC, this->mrmlScene());
  /*vtkMRMLScene* mrmlScene = this->mrmlScene();
  vtkMRMLModelNode* sampleNode = vtkMRMLModelNode::New();
  sampleNode->SetScene(mrmlScene);
  sampleNode->SetName("Sample");
  sampleNode->SetAndObservePolyData(samplePC);
  
  vtkMRMLModelDisplayNode* modelDisplayNode = vtkMRMLModelDisplayNode::New();
  //modelDisplayNode->SetColor(0,1,0); // green
  modelDisplayNode->SetScene(mrmlScene);
  mrmlScene->AddNode(modelDisplayNode);
  sampleNode->SetAndObserveDisplayNodeID(modelDisplayNode->GetID());
  
  modelDisplayNode->SetInputPolyData(sampleNode->GetPolyData());
  mrmlScene->AddNode(sampleNode);*/
  
}

//-----------------------------------------------------------------------------
void qSlicerLoadableSSMBuildingModuleWidget::setMRMLScene(vtkMRMLScene* mrmlScene)
{
  this->qSlicerAbstractModuleWidget::setMRMLScene(mrmlScene);
}

//-----------------------------------------------------------------------------
void qSlicerLoadableSSMBuildingModuleWidget::displayEigenSpectrum()
{
  Q_D(qSlicerLoadableSSMBuildingModuleWidget);
  std::vector<statismo::ScalarType> eigenvalueVector;
  int nbPrincipalComponent;
  double sumEigenvalue = 0;
  // Get the model name selected by the user
  if (d->radioButtonVTK->isChecked()){
    VtkStatisticalModelType* statModel = d->getVtkStatModel();
    VectorType eigenvalue = statModel->GetPCAVarianceVector();
    sumEigenvalue = eigenvalue.sum();
    for (int i=0;i<eigenvalue.rows();i++){
      eigenvalueVector.push_back(eigenvalue[i]);
    }
    nbPrincipalComponent = statModel->GetNumberOfPrincipalComponents();
  }
  if (d->radioButtonITK->isChecked()){
    ItkStatisticalModelType* statModel = d->getItkStatModel();
    itkVectorType itkEigenvalue = statModel->GetPCAVarianceVector(); 
    sumEigenvalue = itkEigenvalue.sum();
    for (int i=0;i<itkEigenvalue.size();i++){
      eigenvalueVector.push_back(itkEigenvalue(i));
    }
    nbPrincipalComponent = statModel->GetNumberOfPrincipalComponents();
  }
  // Switch to a layout (24) that contains a Chart View to initiate the construction of the widget and Chart View Node 
  vtkMRMLScene* mrmlScene = this->mrmlScene();
  mrmlScene->InitTraversal();
  vtkMRMLLayoutNode* sceneLayoutNode = vtkMRMLLayoutNode::SafeDownCast(mrmlScene->GetNextNodeByClass("vtkMRMLLayoutNode") );
  sceneLayoutNode->SetViewArrangement(24);

  // Get the Chart View Node
  vtkMRMLChartViewNode* chartViewNode = vtkMRMLChartViewNode::SafeDownCast(mrmlScene->GetNextNodeByClass("vtkMRMLChartViewNode") );

  // Create an Array Node and add the eigen value
  //VectorType eigenvalue = statModel->GetPCAVarianceVector();
  vtkMRMLDoubleArrayNode* doubleArrayNode = vtkMRMLDoubleArrayNode::SafeDownCast(mrmlScene->AddNode(vtkMRMLDoubleArrayNode::New()));
  vtkDoubleArray* a = doubleArrayNode->GetArray();
  //int nbPrincipalComponent = statModel->GetNumberOfPrincipalComponents();
  std::cout<<"nbPrincipalComponent = "<<nbPrincipalComponent<<std::endl;
  a->SetNumberOfTuples(nbPrincipalComponent);
  double cpt = 0;
  for (int i = 0;i<nbPrincipalComponent;i++){
    cpt = cpt + eigenvalueVector[i]*100/sumEigenvalue;
    a->SetComponent(i, 0, i);
    //a->SetComponent(i, 1, eigenvalueVector[i]);
    a->SetComponent(i, 1, cpt);
    a->SetComponent(i, 2, 0);
  }

  // Create a Chart Node.
  vtkMRMLChartNode* chartNode = vtkMRMLChartNode::New();
  chartNode->AddArray("EigenSpectrum", doubleArrayNode->GetID());
  mrmlScene->AddNode(chartNode);
  
  // Set properties on the Chart
  chartNode->SetProperty("default", "title", "EigenSpectrum");
  chartNode->SetProperty("default", "xAxisLabel", "Principal components");
  chartNode->SetProperty("default", "yAxisLabel", "Eigen value");
  chartNode->SetProperty("default", "type", "Line");
  chartNode->SetProperty("default", "showMarkers", "on");
  
  // Tell the Chart View which Chart to display
  chartViewNode->SetChartNodeID(chartNode->GetID());
}

//-----------------------------------------------------------------------------
void qSlicerLoadableSSMBuildingModuleWidget::onPerformedRegistration()
{
  Q_D(qSlicerLoadableSSMBuildingModuleWidget);
  // we take an arbitrary dataset as the reference, as they have all the same resolution anyway
	 MeshFileReaderType::Pointer refReader = MeshFileReaderType::New();
	 refReader->SetFileName("/home/marine/Documents/SSM/Data/data01/Model_1_1_triangle_dec.vtk");
	 refReader->Update();

	 ItkRepresenterType::Pointer representer = ItkRepresenterType::New();
	 representer->SetReference(refReader->GetOutput());

	 // choose the kernel function that you would like to use for the registration
	 //    statismo::GaussianKernel gk(sigma * sigma); 
  
	 statismo::MultiscaleGaussianKernel mgk(d->sigmaLower->text().toDouble() * d->sigmaLower->text().toDouble(), d->sigmaUpper->text().toDouble() * d->sigmaUpper->text().toDouble(), d->numLevels->text().toInt());
		
  typedef itk::CovarianceFunctionModelBuilder<ItkRepresenterType> ModelBuilderTypeMesh;
	 ModelBuilderTypeMesh::Pointer modelBuilder = ModelBuilderTypeMesh::New();
	 modelBuilder->SetRepresenter(representer);
	 ItkStatisticalModelType::Pointer model = modelBuilder->BuildNewModel(&mgk, d->numPointsForNystrom->text().toInt(), d->numComponents->text().toInt());
	 model->Save("/home/marine/Desktop/vertebraModelFit.h5");
}

//-----------------------------------------------------------------------------
int qSlicerLoadableSSMBuildingModuleWidget::getdir (std::string dir, StringVectorType &files, std::string nameImage, const std::string& extension, std::string &ref, std::string numRef)
{
	itk::Directory::Pointer directory = itk::Directory::New();
	directory->Load(dir.c_str());

	for (unsigned i = 0; i < directory->GetNumberOfFiles(); i++) {
		const char* filename = directory->GetFile(i);
		std::cout<<filename<<std::endl;
		if (std::string(filename).find(extension) != std::string::npos && std::string(filename).find(nameImage) != std::string::npos){
			if (std::string(filename).find(numRef) != std::string::npos){
				ref = dir + "/" +filename;
			}
			files.push_back(dir + "/" +filename);
		}
	}

	std::sort(files.begin(), files.end());
    return 0;
}

//-----------------------------------------------------------------------------
std::vector<PointType> qSlicerLoadableSSMBuildingModuleWidget::readLandmarks(const std::string& filename) {

	std::vector<PointType> ptList;

	std::fstream file ( filename.c_str(), std::fstream::in );
	if (!file) {
			std::cout << "could not read landmark file " << std::endl;
			throw std::runtime_error("could not read landmark file ");
	}
	std::string line;
	while (  std::getline ( file, line))
	{
			if (line.length() > 0 && line[0] == '#')
					continue;

			std::istringstream strstr(line);
			std::string token;
			std::getline(strstr, token, ','); // ignore
			std::getline(strstr, token, ',');
			double pt0 = -atof(token.c_str()); // compensate slicer coordinate system
			std::getline(strstr, token, ',');
			double pt1 = -atof(token.c_str()); // compensate slicer coordinate system
			std::getline(strstr, token, ',');
			double pt2 = atof(token.c_str());
			PointType pt;
			pt[0] = pt0; pt[1] = pt1; pt[2] = pt2;
			ptList.push_back(pt);
	}
	return ptList;
}

//-----------------------------------------------------------------------------
void qSlicerLoadableSSMBuildingModuleWidget::onBuildPartiallyFixedModel()
{ 
  Q_D(qSlicerLoadableSSMBuildingModuleWidget);

  typedef itk::PartiallyFixedModelBuilder<ItkRepresenterType> PartiallyFixedModelBuilderType;
  //typedef itk::PointsLocator<int, 3, double, MeshType::PointsContainer> itkPointsLocatorType;
  typedef itk::PointsLocator<MeshType::PointsContainer> itkPointsLocatorType;
  
  std::string inputModelName="/home/marine/Desktop/vertebraModelFit.h5";
	 std::string dir = "/home/marine/Documents/SSM/Data/Results/Alignment";
	 std::string nameLandmarks = "landmark";
	 std::string nbRef ="01";
	 std::string extension = ".acsv";
	 //float sigma2 = 20;
	 std::string constraintModelName = "/home/marine/Documents/SSM/Data/Results/Registration/vertebraModel";

	 std::string refLandmarksFilename;

	 StringVectorType filenames;
	 getdir(dir, filenames, nameLandmarks, extension, refLandmarksFilename, nbRef);

	 try {
		 // load the model
		 ItkStatisticalModelType::Pointer inputModel = ItkStatisticalModelType::New();
		 inputModel->Load(inputModelName.c_str());

		 int cpt = 0;
		 for (StringVectorType::const_iterator it = filenames.begin(); it != filenames.end(); it++)
		 {
			 PartiallyFixedModelBuilderType::Pointer pfmb = PartiallyFixedModelBuilderType::New();

			 std::vector<PointType> refLandmarks = readLandmarks(refLandmarksFilename);
			 std::vector<PointType> targetLandmarks = readLandmarks(*it);
			 std::cout<<*it<<""<<cpt<<std::endl;

			 // Create an empty list, holding the constraint and add the point
			 ItkStatisticalModelType::PointValueListType constraints;
			 MeshType::Pointer reference = inputModel->GetRepresenter()->GetReference();

			 itkPointsLocatorType::Pointer ptLocator = itkPointsLocatorType::New();
			 ptLocator->SetPoints(reference->GetPoints());
			 ptLocator->Initialize();
			 assert(refLandmarks.size() == targetLandmarks.size());
			 for (unsigned i =0; i < refLandmarks.size(); i++) {
				 // project to closest point on the mesh
				 int closestPointId = ptLocator->FindClosestPoint(refLandmarks[i]);
				 PointType refPoint = (*reference->GetPoints())[closestPointId];
				 std::cout << "refPoint " << refPoint << std::endl;
				 std::cout << "refLm " <<refLandmarks[i] << std::endl;


				 ItkStatisticalModelType::PointValuePairType pointValue(refPoint, targetLandmarks[i]);
				 constraints.push_back(pointValue);
			 }

			 ItkStatisticalModelType::Pointer constraintModel =  pfmb->BuildNewModelFromModel(inputModel, constraints, d->sigma2->text().toFloat());

			 std::ostringstream ss;
			 ss <<constraintModelName<<cpt<<".h5";
			 std::string fullpath = ss.str();
			 constraintModel->Save(fullpath.c_str());
			 cpt++;

			 std::cout << "successfully saved the model to " << constraintModelName << std::endl;
		 }

	 }
	 catch (StatisticalModelException& e) {
		 std::cout << "Exception occurred while building the intensity model" << std::endl;
		 std::cout << e.what() << std::endl;
	 }
}

//-----------------------------------------------------------------------------
  
typedef itk::Image<float, Dimensions> ImageType;
typedef itk::ImageFileReader<ImageType> ImageReaderType;
typedef itk::PenalizingMeanSquaresPointSetToImageMetric<MeshType, ImageType> MetricType;

// As a transform, we use the StatisticalShapeModelTransform, that comes with statismo
typedef itk::StatisticalShapeModelTransform<ItkRepresenterType, double, Dimensions> TransformType;

typedef itk::PointSetToImageRegistrationMethod<MeshType, ImageType> RegistrationFilterType;

typedef  itk::LBFGSOptimizer OptimizerType;

//typedef itk::StatisticalModel<RepresenterTypeMesh> StatisticalModelType;

typedef itk::LinearInterpolateImageFunction<ImageType, double> InterpolatorType;

class IterationStatusObserver : public itk::Command
{
public:
  typedef  IterationStatusObserver   Self;
  typedef  itk::Command             Superclass;
  typedef  itk::SmartPointer<Self>  Pointer;

  itkNewMacro( Self );

  typedef itk::LBFGSOptimizer    OptimizerType;

  typedef const OptimizerType                     *OptimizerPointer;


   void Execute(itk::Object *caller, const itk::EventObject & event)
  {
    Execute( (const itk::Object *)caller, event);
  }

  void Execute(const itk::Object * object, const itk::EventObject & event)
  {
    OptimizerPointer optimizer =
                         dynamic_cast< OptimizerPointer >( object );

    if( ! itk::IterationEvent().CheckEvent( &event ) )
    {
      return;
    }

    std::cout << "Iteration: " << ++m_iter_no ;
    std::cout << "; Value: " << optimizer->GetCachedValue();
    std::cout << "; Current Parameters: " << optimizer->GetCachedCurrentPosition() << std::endl;
  }


protected:
  IterationStatusObserver():
     m_iter_no(0)     {};

  virtual ~IterationStatusObserver(){};

private:
  int m_iter_no;

};

int qSlicerLoadableSSMBuildingModuleWidget::getdir2 (std::string dir, StringVectorType &files, std::string nameImage, const std::string& extension)
{
	itk::Directory::Pointer directory = itk::Directory::New();
	directory->Load(dir.c_str());

	for (unsigned i = 0; i < directory->GetNumberOfFiles(); i++) {
		const char* filename = directory->GetFile(i);

		if (std::string(filename).find(extension) != std::string::npos && std::string(filename).find(nameImage) != std::string::npos){
			files.push_back(dir + "/" +filename);
			//std::cout<<dir + "/" +filename<<std::endl;
		}
	}

	std::sort(files.begin(), files.end());

    return 0;
}


//-----------------------------------------------------------------------------
void qSlicerLoadableSSMBuildingModuleWidget::onShapeModelFitting()
{ 
  Q_D(qSlicerLoadableSSMBuildingModuleWidget);


  std::string dirModel = "/home/marine/Documents/SSM/Data/Results/Registration";
	 std::string modelName = "vertebraModel";
	 std::string extensionModel = ".h5";
	 std::string dirTarget = "/home/marine/Documents/SSM/Data/Results/Alignment";
	 std::string modelTarget = "distance";
	 std::string extensionTarget = ".vtk";
	 //double regularizationParameter = 0.001;
	 std::string outputmeshname = "/home/marine/Documents/SSM/Data/Results/Registration/modelFitted";

	 StringVectorType filenamesModel;
	 getdir2(dirModel, filenamesModel, modelName, extensionModel);

	 StringVectorType filenamesTarget;
	 getdir2(dirTarget, filenamesTarget, modelTarget, extensionTarget);

	 assert(filenamesModel.size() == filenamesTarget.size());

	 for (int i=0; i<filenamesModel.size(); i++)
	 {
		 // load the model
		 ItkStatisticalModelType::Pointer model = ItkStatisticalModelType::New();
		 std::cout<<filenamesModel[i]<<""<<i<<std::endl;
		 model->Load(filenamesModel[i].c_str());
		 MeshType::Pointer fixedPointSet  = model->GetRepresenter()->GetReference();
		 std::cout << "model successfully loaded " << std::endl;

		 // load the image to which we will fit
		 ImageReaderType::Pointer targetReader = ImageReaderType::New();
		 std::cout<<filenamesTarget[i]<<std::endl;
		 targetReader->SetFileName(filenamesTarget[i].c_str());
		 targetReader->Update();
		 ImageType::Pointer movingImage = targetReader->GetOutput();


		 // now we perform the fitting, using the itk registration framework
		 TransformType::Pointer transform = TransformType::New();
		 transform->SetStatisticalModel(model);
		 transform->SetIdentity();

		 // Setting up the fitting
		 OptimizerType::Pointer optimizer = OptimizerType::New();
		 optimizer->SetMaximumNumberOfFunctionEvaluations(1000);
		 optimizer->MinimizeOn();

		 typedef  IterationStatusObserver ObserverType;
		 ObserverType::Pointer observer = ObserverType::New();
		 optimizer->AddObserver( itk::IterationEvent(), observer );

		 MetricType::Pointer metric = MetricType::New();
		 metric->SetRegularizationParameter(d->regularizationParameter->text().toDouble());

		 InterpolatorType::Pointer interpolator = InterpolatorType::New();

		 RegistrationFilterType::Pointer registration = RegistrationFilterType::New();
		 registration->SetInitialTransformParameters(transform->GetParameters());
		 registration->SetInterpolator(interpolator);
		 registration->SetMetric(metric);
		 registration->SetOptimizer(   optimizer);
		 registration->SetTransform(   transform );
		 registration->SetFixedPointSet(  fixedPointSet );
		 registration->SetMovingImage( movingImage);


		 try {
			 std::cout << "starting model fitting" << std::endl;
			 registration->Update();

		 } catch ( itk::ExceptionObject& o ) {
			 std::cout << "caught exception " << o << std::endl;
		 }

		 // We obtain the fitting result by drawing the model instance that belongs to the
		 // optimal tranform parameters (coefficients)
		 MeshType::Pointer mesh = model->DrawSample(transform->GetCoefficients());

		 // Write out the fitting result
		 std::ostringstream ss;
		 ss <<outputmeshname<<i<<extensionTarget;
		 std::string fullpath = ss.str();
		 itk::MeshFileWriter<MeshType>::Pointer writer = itk::MeshFileWriter<MeshType>::New();
		 writer->SetFileName(fullpath.c_str());
		 writer->SetInput(mesh);
		 writer->Update();
	 }
}
