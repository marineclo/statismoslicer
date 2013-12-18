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
#include <itkTransformFileWriter.h>
#include <itkTransformFileReader.h>
//#include <vtkImageToImageFilter.h>
#include <itkVTKImageImport.h>
#include <vtkImageExport.h>
#include <vtkImageData.h>

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
#include "vtkCollection.h"

//Compare PolyData
#include "vtkPointLocator.h"
#include "vtkMath.h"
#include "vtkFloatArray.h"
#include "vtkPointData.h"

//MRML includes
#include <vtkNew.h>
#include "vtkMRMLModelNode.h"
#include "vtkMRMLModelDisplayNode.h"
#include "vtkMRMLScene.h"
#include "vtkMRMLMarkupsFiducialNode.h"
#include "vtkMRMLVolumeNode.h"

#include "qSlicerApplication.h"

const unsigned Dimensions = 3;
typedef itk::Mesh<float, Dimensions  > MeshType;
typedef itk::Point<double, 3> PointType;
typedef itk::PointSet<float, Dimensions > PointSetType;
typedef itk::Image<int, Dimensions> CTImageType;
typedef itk::Image<float, Dimensions> DistanceImageType;
typedef itk::ImageFileReader<CTImageType> CTImageReaderType;
typedef itk::ImageFileWriter<CTImageType> DistanceImageWriterType;
typedef itk::MeshRepresenter<float, Dimensions> RepresenterType;
typedef itk::MeshFileReader<MeshType> MeshReaderType;
typedef itk::MeshFileWriter<MeshType> MeshWriterType;
typedef itk::Image< itk::CovariantVector<float, Dimensions>, Dimensions >   GradientImageType;
//typedef itk::MeanSquaresPointSetToImageMetric<MeshType, DistanceImageType> MetricType;
//typedef itk::PenalizingMeanSquaresPointSetToImageMetric<MeshType, DistanceImageType> MetricType;
typedef itk::PenalizingMeanSquaresPointSetToImageMetric<PointSetType, DistanceImageType> MetricType;
typedef itk::StatisticalShapeModelTransform<RepresenterType, double, Dimensions> StatisticalModelTransformType;
typedef itk::StatisticalModel<RepresenterType> StatisticalModelType;
//typedef itk::PointSetToImageRegistrationMethod<MeshType, DistanceImageType> RegistrationFilterType;
typedef itk::PointSetToImageRegistrationMethod<PointSetType, DistanceImageType> RegistrationFilterType;
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
  typedef const OptimizerType   *OptimizerPointer;

  typedef itk::CompositeTransform<double, 3> CompositeTransformType;
  //typedef const CompositeTransformType   *CompositeTransformPointer;

  void SetParameters(StatisticalModelType::Pointer model, StatisticalModelType::Pointer constraintModel, RigidTransformType::Pointer rigidTransform, vtkMRMLScene* scene, vtkMRMLModelNode* modelNode){
    m_model = model;
    m_constraintModel = constraintModel;
    m_rigidTransform = rigidTransform;
    m_scene = scene;
    m_modelNode = modelNode;
  }

  vtkPolyData* convertMeshToVtk(MeshType::Pointer meshToConvert, vtkPolyData *  m_PolyDataReturn)
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
    //std::cout<<"numPoints = "<<numPoints<<std::endl;

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

    StatisticalModelTransformType::Pointer statModelTransform = StatisticalModelTransformType::New();
    statModelTransform->SetStatisticalModel(m_constraintModel);
    statModelTransform->SetParameters(optimizer->GetCachedCurrentPosition());

    // compose the two transformation
    CompositeTransformType::Pointer transform = CompositeTransformType::New();
    transform->AddTransform(m_rigidTransform);
    transform->AddTransform(statModelTransform);

    //typedef itk::TransformMeshFilter<MeshType, MeshType, CompositeTransformType> TransformMeshFilterType;
    TransformMeshFilterType::Pointer transformMeshFilter = TransformMeshFilterType::New();
    transformMeshFilter->SetInput(m_model->GetRepresenter()->GetReference());
    transformMeshFilter->SetTransform(transform);
    transformMeshFilter->Update();

    vtkPolyData* outputPolyData = vtkPolyData::New();
    outputPolyData=convertMeshToVtk(transformMeshFilter->GetOutput(), outputPolyData);
    vtkSmartPointer< vtkPoints >   point = vtkSmartPointer< vtkPoints >::New();
    point->SetNumberOfPoints(outputPolyData->GetNumberOfPoints());
    for (vtkIdType i=0; i<outputPolyData->GetNumberOfPoints(); i++){
        double pt[3];
        outputPolyData->GetPoint(i,pt);
        pt[0] = -pt[0];
        pt[1] = -pt[1];
        pt[2] = pt[2];
        point->SetPoint(i,pt);
      }
    outputPolyData->SetPoints(point);

   if (m_iter_no>1){
      vtkSmartPointer<vtkCollection> modelDisplayNodes = vtkSmartPointer<vtkCollection>::Take( m_scene->GetNodesByClass("vtkMRMLModelDisplayNode") );
      std::cout<<"Nb Collection= "<<modelDisplayNodes->GetNumberOfItems ()<<std::endl;
      vtkMRMLModelDisplayNode* modelViewNode = vtkMRMLModelDisplayNode::SafeDownCast( modelDisplayNodes->GetItemAsObject (modelDisplayNodes->GetNumberOfItems ()-1));
      modelViewNode->VisibilityOff();
    }


    //Add to the scene
    vtkNew<vtkMRMLModelDisplayNode> modelDisplayNode;
    m_scene->AddNode(modelDisplayNode.GetPointer());

    QString outputModelName = "OutputModelSeg" +QString::number(m_iter_no);

   /* m_modelNode->SetAndObservePolyData(outputPolyData);
    m_modelNode->SetAndObserveDisplayNodeID(modelDisplayNode->GetID());
    m_modelNode->SetName(outputModelName.toStdString().c_str());
    m_modelNode->Modified();*/

    vtkNew<vtkMRMLModelNode> meanNode;
    meanNode->SetAndObservePolyData(outputPolyData);
    meanNode->SetAndObserveDisplayNodeID(modelDisplayNode->GetID());
    meanNode->SetName(outputModelName.toStdString().c_str());
    m_scene->AddNode(meanNode.GetPointer());

    qSlicerApplication * app = qSlicerApplication::application();
    app->processEvents();


    outputPolyData->Delete();
  }


protected:
  IterationStatusObserver():
     m_iter_no(0)     {};

  virtual ~IterationStatusObserver(){};

  StatisticalModelType::Pointer m_model;
  StatisticalModelType::Pointer m_constraintModel;
  RigidTransformType::Pointer m_rigidTransform;
  vtkMRMLScene* m_scene;
  vtkMRMLModelNode* m_modelNode;

private:
  int m_iter_no;

};


//-----------------------------------------------------------------------------
/// \ingroup Slicer_QtModules_ExtensionTemplate
class qSlicerLandmarkSegmentationModuleWidgetPrivate: public Ui_qSlicerLandmarkSegmentationModuleWidget
{
public:
  qSlicerLandmarkSegmentationModuleWidgetPrivate();
  
  std::vector<PointType > readLandmarks(vtkMRMLMarkupsFiducialNode* markupsNode);
  StatisticalModelType::Pointer computePartiallyFixedModel(const RigidTransformType* rigidTransform, const StatisticalModelType* statisticalModel,
const  std::vector<PointType >& modelLandmarks,const  std::vector<PointType >& targetLandmarks, double variance);

  vtkPolyData* convertMeshToVtk(MeshType::Pointer meshToConvert, vtkPolyData *  m_PolyDataReturn);
  
  void ConnectVTKToITK(vtkImageExport* in, itk::VTKImageImport<CTImageType>* out);
};

//-----------------------------------------------------------------------------
// qSlicerLandmarkSegmentationModuleWidgetPrivate methods

//-----------------------------------------------------------------------------
qSlicerLandmarkSegmentationModuleWidgetPrivate::qSlicerLandmarkSegmentationModuleWidgetPrivate()
{
}

/**
* read landmarks from Slicer list and return them as a list.
*
* The format is: label,x,y,z
*
* @param Markups Node
* @returns A list of itk points
*/
std::vector<PointType > qSlicerLandmarkSegmentationModuleWidgetPrivate::readLandmarks(vtkMRMLMarkupsFiducialNode* markupsNode) {

  std::vector<PointType> ptList;

  for (int i=0;i<markupsNode->GetNumberOfFiducials();i++){
      double pos[3];
      markupsNode->GetNthFiducialPosition(i,pos);
      PointType pt;
      pt[0] = -pos[0];
      pt[1] = -pos[1];
      pt[2] = pos[2];
      ptList.push_back(pt);
      //std::cout<<"pt[0]= "<<pt[0]<<" pt[1]= "<<pt[1]<<" pt[2]= "<<pt[2]<<std::endl;
    }

  return ptList;
}

// Returns a new model, that is restricted to go through the points specified in targetLandmarks..
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
  //std::cout<<"numPoints = "<<numPoints<<std::endl;

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

void qSlicerLandmarkSegmentationModuleWidgetPrivate::ConnectVTKToITK(vtkImageExport* in, itk::VTKImageImport<CTImageType>* out){
  
  out->SetUpdateInformationCallback(in->GetUpdateInformationCallback());
  out->SetPipelineModifiedCallback(in->GetPipelineModifiedCallback());
  out->SetWholeExtentCallback(in->GetWholeExtentCallback());
  out->SetSpacingCallback(in->GetSpacingCallback());
  out->SetOriginCallback(in->GetOriginCallback());
  out->SetScalarTypeCallback(in->GetScalarTypeCallback());
  out->SetNumberOfComponentsCallback(in->GetNumberOfComponentsCallback());
  out->SetPropagateUpdateExtentCallback(in->GetPropagateUpdateExtentCallback());
  out->SetUpdateDataCallback(in->GetUpdateDataCallback());
  out->SetDataExtentCallback(in->GetDataExtentCallback());
  out->SetBufferPointerCallback(in->GetBufferPointerCallback());
  out->SetCallbackUserData(in->GetCallbackUserData());

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

void qSlicerLandmarkSegmentationModuleWidget::setModel(){
  Q_D(qSlicerLandmarkSegmentationModuleWidget);
  QString inputFile = QFileDialog::getOpenFileName(this, "Select model", QString());
  d->modelName->setText(inputFile);
}

/*void qSlicerLandmarkSegmentationModuleWidget::setCTscan(){
  Q_D(qSlicerLandmarkSegmentationModuleWidget);
  QString inputFile = QFileDialog::getOpenFileName(this, "Select CT scan", QString());
  d->targetName->setText(inputFile);
}*/


void qSlicerLandmarkSegmentationModuleWidget::apply(){
  Q_D(qSlicerLandmarkSegmentationModuleWidget);
  
  vtkMRMLVolumeNode* volumeNode = vtkMRMLVolumeNode::SafeDownCast(d->qMRMLVolumeNodeComboBox->currentNode());
  vtkImageData* vtkImage = volumeNode->GetImageData();
  
  //vtkImageDataInput->SetScalarTypeToFloat();
  double pt[3];
  double pt1[3];
  volumeNode->GetSpacing(pt);
  volumeNode->GetOrigin(pt1);
  vtkImageData* vtkImageDataInput = vtkImageData::New();
  vtkImageDataInput->DeepCopy(vtkImage);
  vtkImageDataInput->SetSpacing(pt);
  vtkImageDataInput->SetOrigin(pt1);
  
  typedef itk::VTKImageImport<CTImageType> InputImageImportType;
  vtkImageExport* inputImageExporter = vtkImageExport::New();
  inputImageExporter->SetInput(vtkImageDataInput);
  InputImageImportType::Pointer inputImageImporter = InputImageImportType::New();

  d->ConnectVTKToITK(inputImageExporter, inputImageImporter);
  CTImageType* ctImage = const_cast<CTImageType*>(inputImageImporter->GetOutput());
  ctImage->Update();
  
 /* DistanceImageWriterType::Pointer targetWriter = DistanceImageWriterType::New();
  targetWriter->SetFileName("/tmp/volume.nrrd");
  targetWriter->SetInput(ctImage);
  targetWriter->Update();*/
  
  
  // load the image to which we will fit
  /*CTImageReaderType::Pointer targetReader = CTImageReaderType::New();
  targetReader->SetFileName(d->targetName->text().toStdString().c_str());
  targetReader->Update();
  CTImageType::Pointer ctImage = targetReader->GetOutput();*/

  // We compute a binary threshold of input image to get a rough segmentation of the bony structure.
  // Then we compute a distance transform of the segmentation, which we then use for the fitting
  BinaryThresholdImageFilterType::Pointer thresholdFilter = BinaryThresholdImageFilterType::New();
  thresholdFilter->SetInput(ctImage);
  thresholdFilter->SetLowerThreshold( d->threshold->text().toFloat() );
  CTImageType::Pointer thresholdImage = thresholdFilter->GetOutput();
  
  DistanceMapImageFilterType::Pointer dm = DistanceMapImageFilterType::New();
  dm->SetInput(thresholdImage);
  dm->Update();
  DistanceImageType::Pointer distanceImage = dm->GetOutput();

  // read the landmarks
  vtkMRMLMarkupsFiducialNode* movingMarkupsNode =  vtkMRMLMarkupsFiducialNode::SafeDownCast(d->qMRMLMarkupsNodeComboBox->currentNode());
  vtkMRMLMarkupsFiducialNode* fixedMarkupsNode =  vtkMRMLMarkupsFiducialNode::SafeDownCast(d->qMRMLRefMarkupsNodeComboBox->currentNode());
  std::vector<PointType> fixedLandmarks = d->readLandmarks(fixedMarkupsNode);
  std::vector<PointType> movingLandmarks = d->readLandmarks(movingMarkupsNode);

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

  vtkMRMLModelNode* modelNode = vtkMRMLModelNode::SafeDownCast(d->qMRMLModelNodeComboBox->currentNode());
  observer->SetParameters(model, constraintModel, rigidTransform, this->mrmlScene(),modelNode);

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

  // we create the fixedPointSet of the registration from the reference mesh of our model.
  // As we are fitting to the 0 level set of a distance image, we set the value of the pointdata to 0.
  PointSetType::Pointer fixedPointSet = PointSetType::New();
  fixedPointSet->SetPoints(model->GetRepresenter()->GetReference()->GetPoints());
  PointSetType::PointDataContainer::Pointer points = PointSetType::PointDataContainer::New();
  points->Reserve(model->GetRepresenter()->GetReference()->GetNumberOfPoints());
  for (PointSetType::PointDataContainer::Iterator it = points->Begin(); it != points->End(); ++it) {
      it->Value() = 0;
  }
  fixedPointSet->SetPointData(points);
  registration->SetFixedPointSet(  fixedPointSet);
  //registration->SetFixedPointSet(  model->GetRepresenter()->GetReference() );
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
  
  vtkPolyData* outputPolyData = vtkPolyData::New();
  outputPolyData=d->convertMeshToVtk(transformMeshFilter->GetOutput(), outputPolyData);
  
  //Add to the scene
  vtkNew<vtkMRMLModelDisplayNode> modelDisplayNode;
  this->mrmlScene()->AddNode(modelDisplayNode.GetPointer());

  vtkNew<vtkMRMLModelNode> meanNode;
  meanNode->SetAndObservePolyData(outputPolyData);
  meanNode->SetAndObserveDisplayNodeID(modelDisplayNode->GetID());
  meanNode->SetName("OutputModelSegFinal");
  this->mrmlScene()->AddNode(meanNode.GetPointer());
  
  outputPolyData->Delete();
      
}

//-----------------------------------------------------------------------------
void qSlicerLandmarkSegmentationModuleWidget::comparePolyData()
{
  Q_D(qSlicerLandmarkSegmentationModuleWidget);

  vtkMRMLModelNode* model1Node = vtkMRMLModelNode::SafeDownCast(d->qMRMLModel1NodeComboBox->currentNode());
  vtkMRMLModelNode* model2Node = vtkMRMLModelNode::SafeDownCast(d->qMRMLModel2NodeComboBox->currentNode());
  vtkSmartPointer< vtkPolyData > polyData1 = model1Node->GetPolyData();
  vtkSmartPointer< vtkPolyData > polyData2 = model2Node->GetPolyData();

  vtkSmartPointer<vtkPointLocator> pointLocator2 = vtkSmartPointer<vtkPointLocator>::New();
  pointLocator2->SetDataSet(polyData2);

  vtkSmartPointer<vtkFloatArray> scalars = vtkSmartPointer<vtkFloatArray>::New();
  scalars->SetNumberOfComponents(1);
  scalars->Allocate(polyData1->GetNumberOfPoints(),1000);

  for (vtkIdType i=0; i<polyData1->GetNumberOfPoints(); i++){
      double Xi[3];
      polyData1->GetPoint(i,Xi);
      double Yid[3];
      int ptLocatorId = pointLocator2->FindClosestPoint(Xi);
      polyData2->GetPoint(ptLocatorId,Yid);
      double distance = sqrt(vtkMath::Distance2BetweenPoints(Xi,Yid));
      scalars->InsertTuple1(i,distance);
    }

  vtkSmartPointer<vtkPointLocator> pointLocator1 = vtkSmartPointer<vtkPointLocator>::New();
  pointLocator1->SetDataSet(polyData1);

  for (vtkIdType i=0; i<polyData2->GetNumberOfPoints(); i++){
      double Yi[3];
      polyData2->GetPoint(i,Yi);
      double Xid[3];
      vtkIdType ptLocatorId = pointLocator1->FindClosestPoint(Yi);
      polyData1->GetPoint(ptLocatorId,Xid);
      double distance = sqrt(vtkMath::Distance2BetweenPoints(Yi,Xid));
      double xdistance = scalars->GetTuple1(ptLocatorId);
      if (distance>xdistance){
          scalars->SetTuple1(ptLocatorId,distance);
        }
    }


  vtkSmartPointer< vtkPolyData > polyData3 = vtkSmartPointer< vtkPolyData >::New();
  polyData3->DeepCopy(polyData1);
  polyData3->GetPointData()->SetScalars(scalars);
  scalars->SetName("scalarsDis");
  double rangeScalars[2];
  scalars->GetRange(rangeScalars);
  std::cout<<"min= "<<rangeScalars[0]<<" max = "<<rangeScalars[1]<<std::endl;

  // Add polydata to the scene
  vtkNew<vtkMRMLModelDisplayNode> sampleDisplayNode;
  sampleDisplayNode->SetScalarVisibility(true);
  sampleDisplayNode->SetActiveScalarName(scalars->GetName());
  sampleDisplayNode->SetAutoScalarRange(true);
  sampleDisplayNode->SetAndObserveColorNodeID("vtkMRMLColorTableNodeWarm1");
  this->mrmlScene()->AddNode(sampleDisplayNode.GetPointer());

  vtkNew<vtkMRMLModelNode> sampleNode;
  sampleNode->SetAndObservePolyData(polyData3);
  sampleNode->SetAndObserveDisplayNodeID(sampleDisplayNode->GetID());

  sampleNode->SetName(d->outputPolyDataName->text().toStdString().c_str());
  this->mrmlScene()->AddNode(sampleNode.GetPointer());

}
