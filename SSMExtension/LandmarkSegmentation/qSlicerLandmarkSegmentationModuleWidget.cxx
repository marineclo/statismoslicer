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
#include "qSlicerLandmarkSegmentationModuleWidget.h"
#include "ui_qSlicerLandmarkSegmentationModuleWidget.h"

// Statismo includes
#include "itkMeshRepresenter.h"
#include "statismo_ITK/itkStatisticalModel.h"
#include "statismo_ITK/itkStatisticalShapeModelTransform.h"
#include "statismo_ITK/itkPosteriorModelBuilder.h"

// ITK includes
#include "itkCommand.h"
#include "itkMesh.h"
#include "itkMeshFileReader.h"
#include "itkMeshFileWriter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
//#include "itkMeanSquaresPointSetToImageMetric.h"
#include "itkPenalizingMeanSquaresPointSetToImageMetric.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkLBFGSOptimizer.h"
#include "itkPointSetToImageRegistrationMethod.h"
#include "itkTransformMeshFilter.h"
#include "itkCompositeTransform.h"
#include "itkRigid3DTransform.h"
#include "itkLandmarkBasedTransformInitializer.h"
#include "itkCastImageFilter.h"
#include "itkCannyEdgeDetectionImageFilter.h"
#include "itkSignedDanielssonDistanceMapImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkPointsLocator.h"

// convert itk Mesh to vtk PolyData
#include "vtkPolyData.h"
#include "vtkPolyDataWriter.h"
#include "vtkPoints.h"
#include "vtkCellArray.h"
#include "itkDefaultDynamicMeshTraits.h"
#include "itkTriangleCell.h"
#include "itkPoint.h"
#include "itkObject.h"
#include <vtkSmartPointer.h>

//MRML includes
#include <vtkNew.h>
#include "vtkMRMLModelNode.h"
#include "vtkMRMLModelDisplayNode.h"
#include "vtkMRMLScene.h"

const unsigned Dimensions = 3;
typedef itk::Mesh<float, Dimensions  > MeshType;
typedef itk::Point<double, 3> PointType;
typedef itk::Image<float, Dimensions> CTImageType;
typedef itk::Image<float, Dimensions> DistanceImageType;
typedef itk::ImageFileReader<CTImageType> CTImageReaderType;
typedef itk::ImageFileWriter<DistanceImageType> DistanceImageWriterType;
typedef itk::MeshRepresenter<float, Dimensions> RepresenterType;
typedef itk::MeshFileReader<MeshType> MeshReaderType;
typedef itk::Image< itk::CovariantVector<float, Dimensions>, Dimensions >   GradientImageType;
//typedef itk::MeanSquaresPointSetToImageMetric<MeshType, DistanceImageType> MetricType;
typedef itk::PenalizingMeanSquaresPointSetToImageMetric<MeshType, DistanceImageType> MetricType;
typedef itk::StatisticalShapeModelTransform<RepresenterType, double, Dimensions> StatisticalModelTransformType;
typedef itk::StatisticalModel<RepresenterType> StatisticalModelType;
typedef itk::PointSetToImageRegistrationMethod<MeshType, DistanceImageType> RegistrationFilterType;
typedef itk::LBFGSOptimizer OptimizerType;
typedef itk::LinearInterpolateImageFunction<DistanceImageType, double> InterpolatorType;
typedef itk::BinaryThresholdImageFilter <CTImageType, CTImageType>  BinaryThresholdImageFilterType;
typedef itk::SignedDanielssonDistanceMapImageFilter<CTImageType, DistanceImageType> DistanceMapImageFilterType;
typedef itk::VersorRigid3DTransform<double> RigidTransformType;
typedef itk::LandmarkBasedTransformInitializer<RigidTransformType, DistanceImageType, DistanceImageType> LandmarkTransformInitializerType;
typedef itk::CompositeTransform<double, 3> CompositeTransformType;
typedef itk::TransformMeshFilter<MeshType, MeshType, CompositeTransformType> TransformMeshFilterType;
typedef itk::PosteriorModelBuilder<RepresenterType> PartiallyFixedModelBuilderType;

#if (ITK_VERSION_MAJOR == 4 && ITK_VERSION_MINOR >= 4)
typedef itk::PointsLocator< MeshType::PointsContainer > PointsLocatorType;
#else
typedef itk::PointsLocator<int, 3, double, MeshType::PointsContainer > PointsLocatorType;
#endif


//
// This class is used to track the progress of the optimization
// (its method Execute is called in each iteration of the optimization)
//
class IterationStatusObserver : public itk::Command
{
public:
  typedef  IterationStatusObserver   Self;
  typedef  itk::Command             Superclass;
  typedef  itk::SmartPointer<Self>  Pointer;

  itkNewMacro( Self );

  typedef itk::LBFGSOptimizer    OptimizerType;
  //typedef itk::GradientDescentOptimizer OptimizerType;

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


//-----------------------------------------------------------------------------
/// \ingroup Slicer_QtModules_ExtensionTemplate
class qSlicerLandmarkSegmentationModuleWidgetPrivate: public Ui_qSlicerLandmarkSegmentationModuleWidget
{
public:
  qSlicerLandmarkSegmentationModuleWidgetPrivate();
  
  std::vector<PointType > readLandmarks(const std::string& filename);
  StatisticalModelType::Pointer computePartiallyFixedModel(const RigidTransformType* rigidTransform, const StatisticalModelType* statisticalModel,
const  std::vector<PointType >& modelLandmarks,const  std::vector<PointType >& targetLandmarks, double variance);

  vtkPolyData* convertMeshToVtk(MeshType::Pointer meshToConvert, vtkPolyData *  m_PolyDataReturn);
};

//-----------------------------------------------------------------------------
// qSlicerLandmarkSegmentationModuleWidgetPrivate methods

//-----------------------------------------------------------------------------
qSlicerLandmarkSegmentationModuleWidgetPrivate::qSlicerLandmarkSegmentationModuleWidgetPrivate()
{
}

/**
* read landmarks from the given file in slicer fcsv formant and return them as a list.
*
* The format is: label,x,y,z
*
* @param filename the filename
* @returns A list of itk points
*/
std::vector<PointType > qSlicerLandmarkSegmentationModuleWidgetPrivate::readLandmarks(const std::string& filename) {

  std::vector<PointType> ptList;

  std::fstream file ( filename.c_str() );
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
    std::getline(strstr, token, ','); // ignore the label
    std::getline(strstr, token, ','); // get the x coord
    double pt0 = -atof(token.c_str());
    std::getline(strstr, token, ','); // get the y coord
    double pt1 = -atof(token.c_str());
    std::getline(strstr, token, ','); // get the z coord
    double pt2 = atof(token.c_str());
    PointType pt;
    pt[0] = pt0; pt[1] = pt1; pt[2] = pt2;
    ptList.push_back(pt);
  }
  return ptList;
}

// Returns a new model, that is restricted to go through the proints specified in targetLandmarks..
//
StatisticalModelType::Pointer qSlicerLandmarkSegmentationModuleWidgetPrivate::computePartiallyFixedModel(const RigidTransformType* rigidTransform, const StatisticalModelType* statisticalModel,
const  std::vector<PointType >& modelLandmarks,const  std::vector<PointType >& targetLandmarks, double variance)
{

  // invert the transformand back transform the landmarks
  RigidTransformType::Pointer rinv = RigidTransformType::New();
  rigidTransform->GetInverse(rinv);

  StatisticalModelType::PointValueListType constraints;

  // We need to make sure the the points in fixed landmarks are real vertex points of the model reference.
  MeshType::Pointer reference = statisticalModel->GetRepresenter()->GetReference();
  PointsLocatorType::Pointer ptLocator = PointsLocatorType::New();
  ptLocator->SetPoints(reference->GetPoints());
  ptLocator->Initialize();

  assert(modelLandmarks.size() == targetLandmarks.size());
  for (unsigned i = 0; i < targetLandmarks.size(); i++) {

	int closestPointId = ptLocator->FindClosestPoint(modelLandmarks[i]);
    PointType refPoint = (*reference->GetPoints())[closestPointId];

    // compensate for the rigid transformation that was applied to the model
    PointType targetLmAtModelPos = rinv->TransformPoint(targetLandmarks[i]);
    StatisticalModelType::PointValuePairType pointValue(refPoint ,targetLmAtModelPos);
    constraints.push_back(pointValue);

  }

  PartiallyFixedModelBuilderType::Pointer partiallyFixedModelBuilder = PartiallyFixedModelBuilderType::New();
  StatisticalModelType::Pointer partiallyFixedModel = partiallyFixedModelBuilder->BuildNewModelFromModel(statisticalModel,constraints, variance, false);

  return partiallyFixedModel;
}


vtkPolyData* qSlicerLandmarkSegmentationModuleWidgetPrivate::convertMeshToVtk(MeshType::Pointer meshToConvert, vtkPolyData *  m_PolyDataReturn)
{
  typedef MeshType::MeshTraits                      TriangleMeshTraits;
  typedef MeshType::PointType                       PointType;
  typedef MeshType::PointsContainer                 InputPointsContainer;
  typedef InputPointsContainer::Pointer             InputPointsContainerPointer;
  typedef InputPointsContainer::Iterator            InputPointsContainerIterator;
  typedef MeshType::CellType                        CellType;

  typedef MeshType::CellsContainerPointer           CellsContainerPointer;
  typedef MeshType::CellsContainerIterator          CellsContainerIterator;

  //vtkPoints  *   m_Points = vtkPoints::New();
  vtkSmartPointer< vtkPoints >   m_Points = vtkSmartPointer< vtkPoints >::New();
  vtkSmartPointer< vtkPolyData >   m_PolyData = vtkSmartPointer< vtkPolyData >::New();

  vtkCellArray * m_Polys = vtkCellArray::New();

  int numPoints =  meshToConvert->GetNumberOfPoints();
  std::cout<<"numPoints = "<<numPoints<<std::endl;

  InputPointsContainerPointer      myPoints = meshToConvert->GetPoints();
  InputPointsContainerIterator     points = myPoints->Begin();
  PointType point;

  if (numPoints == 0)
  {
    printf( "Aborting: No Points in GRID\n");
  }

  m_Points->SetNumberOfPoints(numPoints);
  
  int idx=0;
  float vpoint[3];
  while( points != myPoints->End() )
  {
    point = points.Value();
    vpoint[0]= point[0];
    vpoint[1]= point[1];
    vpoint[2]= point[2];
    m_Points->SetPoint(idx++,vpoint);
    points++;
  }

  m_PolyData->SetPoints(m_Points);

  //m_Points->Delete();

  CellsContainerPointer cells = meshToConvert->GetCells();
  CellsContainerIterator cellIt = cells->Begin();
  vtkIdType pts[3];
  while ( cellIt != cells->End() )
  {
    CellType *nextCell = cellIt->Value();
    CellType::PointIdIterator pointIt = nextCell->PointIdsBegin();
    int i;

    switch (nextCell->GetType())
    {
      case CellType::VERTEX_CELL:
      case CellType::LINE_CELL:
      case CellType::POLYGON_CELL:
        break;
      case CellType::TRIANGLE_CELL:
        i=0;
        while (pointIt != nextCell->PointIdsEnd() )
        {
          pts[i++] = *pointIt++;
        }
        m_Polys->InsertNextCell(3,pts);
        break;
      default:
        printf("something \n");
    }
    cellIt++;

  }

  m_PolyData->SetPolys(m_Polys);
  m_Polys->Delete();
 
  m_PolyDataReturn->ShallowCopy(m_PolyData);

  return m_PolyDataReturn;
}

//-----------------------------------------------------------------------------
// qSlicerLandmarkSegmentationModuleWidget methods

//-----------------------------------------------------------------------------
qSlicerLandmarkSegmentationModuleWidget::qSlicerLandmarkSegmentationModuleWidget(QWidget* _parent)
  : Superclass( _parent )
  , d_ptr( new qSlicerLandmarkSegmentationModuleWidgetPrivate )
{
}

//-----------------------------------------------------------------------------
qSlicerLandmarkSegmentationModuleWidget::~qSlicerLandmarkSegmentationModuleWidget()
{
}

//-----------------------------------------------------------------------------
void qSlicerLandmarkSegmentationModuleWidget::setup()
{
  Q_D(qSlicerLandmarkSegmentationModuleWidget);
  d->setupUi(this);
  //this->Superclass::setup();
}

//-----------------------------------------------------------------------------
void qSlicerLandmarkSegmentationModuleWidget::setMRMLScene(vtkMRMLScene* mrmlScene)
{
  this->qSlicerAbstractModuleWidget::setMRMLScene(mrmlScene);
}

void qSlicerLandmarkSegmentationModuleWidget::setRefFid(){

  Q_D(qSlicerLandmarkSegmentationModuleWidget);
  QString inputFile = QFileDialog::getOpenFileName(this, "Select Reference Fiducial", QString());
  d->fixedLandmarksName->setText(inputFile);
}

void qSlicerLandmarkSegmentationModuleWidget::setModel(){
  Q_D(qSlicerLandmarkSegmentationModuleWidget);
  QString inputFile = QFileDialog::getOpenFileName(this, "Select model", QString());
  d->modelName->setText(inputFile);
}

void qSlicerLandmarkSegmentationModuleWidget::setInputFid(){
  Q_D(qSlicerLandmarkSegmentationModuleWidget);
  QString inputFile = QFileDialog::getOpenFileName(this, "Select Input Fiducial", QString());
  d->movingLandmarksName->setText(inputFile);
}

void qSlicerLandmarkSegmentationModuleWidget::setCTscan(){
  Q_D(qSlicerLandmarkSegmentationModuleWidget);
  QString inputFile = QFileDialog::getOpenFileName(this, "Select CT scan", QString());
  d->targetName->setText(inputFile);
}


void qSlicerLandmarkSegmentationModuleWidget::apply(){
  Q_D(qSlicerLandmarkSegmentationModuleWidget);

  // load the image to which we will fit
  CTImageReaderType::Pointer targetReader = CTImageReaderType::New();
  targetReader->SetFileName(d->targetName->text().toStdString().c_str());
  targetReader->Update();
  CTImageType::Pointer ctImage = targetReader->GetOutput();

  // We compute a binary threshold of input image to get a rough segmentation of the bony structure.
  // Then we compute a distance transform of the segmentation, which we then use for the fitting
  BinaryThresholdImageFilterType::Pointer thresholdFilter = BinaryThresholdImageFilterType::New();
  thresholdFilter->SetInput(ctImage);
  thresholdFilter->SetLowerThreshold( d->threshold->text().toFloat() );
  CTImageType::Pointer thresholdImage = thresholdFilter->GetOutput();
  
  DistanceImageWriterType::Pointer thresholdWriter = DistanceImageWriterType::New();
  thresholdWriter->SetInput(thresholdFilter->GetOutput());
  //thresholdWriter->SetFileName("/Volumes/MYPASSPORT/TestSegmentation/data136TestSeg-label.nrrd");
  thresholdWriter->SetFileName("/tmp/data136TestSeg-label.nrrd");
  thresholdWriter->Update();
        
  /*CTImageReaderType::Pointer thresholdReader = CTImageReaderType::New();
  thresholdReader->SetFileName(thresholdImageName);
  thresholdReader->Update();
  CTImageType::Pointer thresholdImage = thresholdReader->GetOutput();*/
        
  DistanceMapImageFilterType::Pointer dm = DistanceMapImageFilterType::New();
  dm->SetInput(thresholdImage);
  dm->Update();
  DistanceImageType::Pointer distanceImage = dm->GetOutput();

  // read the landmarks
  std::vector<PointType> fixedLandmarks = d->readLandmarks(d->fixedLandmarksName->text().toStdString().c_str());
  std::vector<PointType> movingLandmarks = d->readLandmarks(d->movingLandmarksName->text().toStdString().c_str());

  // initialize the rigid transform
  RigidTransformType::Pointer rigidTransform = RigidTransformType::New();
  LandmarkTransformInitializerType::Pointer initializer = LandmarkTransformInitializerType::New();
  initializer->SetFixedLandmarks(fixedLandmarks);
  initializer->SetMovingLandmarks(movingLandmarks);
  initializer->SetTransform(rigidTransform);
  initializer->InitializeTransform();

  // load the model create a shape model transform with it
  StatisticalModelType::Pointer model = StatisticalModelType::New();
  model->Load(d->modelName->text().toStdString().c_str());

  StatisticalModelType::Pointer constraintModel = d->computePartiallyFixedModel(rigidTransform, model, fixedLandmarks, movingLandmarks, d->lmVariance->text().toDouble());

  StatisticalModelTransformType::Pointer statModelTransform = StatisticalModelTransformType::New();
  statModelTransform->SetStatisticalModel(constraintModel);
  statModelTransform->SetIdentity();

  // compose the two transformation
  CompositeTransformType::Pointer transform = CompositeTransformType::New();
  transform->AddTransform(rigidTransform);
  transform->AddTransform(statModelTransform);
  transform->SetOnlyMostRecentTransformToOptimizeOn(); //  only optimize the shape parameters, not the rigid transform parameters
  //transform->SetAllTransformsToOptimizeOn(); // optimize shape and pose parameters

  // Some configuration parameters
  short maxNumberOfIterations = 10000; // the maximum number of iterations to use in the optimization
  double translationScale = 1; // dynamic range of translations
  double rotationScale = 0.1; // dynamic range of rotations
  double smScale = 3; // dynamic range of statistical model parameters

  // Setting up the fitting
  OptimizerType::Pointer optimizer = OptimizerType::New();
  optimizer->SetMaximumNumberOfFunctionEvaluations(maxNumberOfIterations);
  optimizer->MinimizeOn();

  unsigned numStatmodelParameters = statModelTransform->GetNumberOfParameters();
  unsigned totalNumParameters =  rigidTransform->GetNumberOfParameters() + numStatmodelParameters;
  
  // set the scales of the optimizer, to compensate for potentially different scales of translation, rotation and shape parameters
  OptimizerType::ScalesType scales( totalNumParameters );
  for (unsigned i = 0; i < numStatmodelParameters; i++) {
    scales[i] = 1.0 / (smScale);
  }
  for (unsigned i = numStatmodelParameters; i < numStatmodelParameters + 3; i++) {
    scales[i] =  1.0 / (rotationScale);
  }
  for (unsigned i = numStatmodelParameters + 3; i < statModelTransform->GetNumberOfParameters() + 6; i++) {
     scales[i] = 1.0 / (translationScale);
  }
  optimizer->SetScales(scales);

  // set up the observer to keep track of the progress
  typedef  IterationStatusObserver ObserverType;
  ObserverType::Pointer observer = ObserverType::New();
  optimizer->AddObserver( itk::IterationEvent(), observer );

  // set up the metric and interpolators
  MetricType::Pointer metric = MetricType::New();
  metric->SetRegularizationParameter(0.01);
  InterpolatorType::Pointer interpolator = InterpolatorType::New();

  // connect all the components
  RegistrationFilterType::Pointer registration = RegistrationFilterType::New();
  registration->SetInitialTransformParameters(transform->GetParameters());
  registration->SetInterpolator(interpolator);
  registration->SetMetric(metric);
  registration->SetOptimizer(   optimizer);
  registration->SetTransform(   transform );


  // the input to the registration will be the reference of the statistical model and the
  // distance map we computed above.
  registration->SetFixedPointSet(  model->GetRepresenter()->GetReference() );
  registration->SetMovingImage(distanceImage);

  try {
	std::cout << "starting model fitting" << std::endl;
    registration->Update();

  } catch ( itk::ExceptionObject& o ) {
    std::cout << "caught exception " << o << std::endl;
  }

  TransformMeshFilterType::Pointer transformMeshFilter = TransformMeshFilterType::New();
  transformMeshFilter->SetInput(model->GetRepresenter()->GetReference());
  transformMeshFilter->SetTransform(transform);
  transformMeshFilter->Update();

  MeshType::Pointer meshOut = transformMeshFilter->GetOutput();
  int numPoints =  meshOut->GetNumberOfPoints();
  std::cout<<"numPoints1 = "<<numPoints<<std::endl;
  
  vtkPolyData* outputPolyData = vtkPolyData::New();
  outputPolyData=d->convertMeshToVtk(transformMeshFilter->GetOutput(), outputPolyData);
  
  //Add to the scene
  vtkNew<vtkMRMLModelDisplayNode> modelDisplayNode;
  this->mrmlScene()->AddNode(modelDisplayNode.GetPointer());

  vtkNew<vtkMRMLModelNode> meanNode;
  meanNode->SetAndObservePolyData(outputPolyData);
  meanNode->SetAndObserveDisplayNodeID(modelDisplayNode->GetID());
  meanNode->SetName("OutputModelSeg");
  this->mrmlScene()->AddNode(meanNode.GetPointer());
  
  // Write out the fitting result
  itk::MeshFileWriter<MeshType>::Pointer writer = itk::MeshFileWriter<MeshType>::New();
  writer->SetFileName("/tmp/mesh136result.vtk");
  writer->SetInput(transformMeshFilter->GetOutput());
  writer->Update();
  
  outputPolyData->Delete();
        
}
