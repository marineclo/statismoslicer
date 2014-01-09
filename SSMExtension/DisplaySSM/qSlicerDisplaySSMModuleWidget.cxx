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
#include "qSlicerDisplaySSMModuleWidget.h"
#include "ui_qSlicerDisplaySSMModuleWidget.h"

// Statismo includes
#include "itkMeshRepresenter.h"
#include "itkMesh.h"
#include "statismo_ITK/itkStatisticalModel.h"

//Add to convert itkMesh to vtkPolyData
#include "vtkPolyData.h"
#include "vtkPolyDataWriter.h"
#include "vtkPoints.h"
#include "vtkCellArray.h"
#include "itkDefaultDynamicMeshTraits.h"
#include "itkTriangleCell.h"
#include "itkPoint.h"
#include "itkObject.h"
#include <vtkSmartPointer.h>

//Add to convert vtkPolyData to vtkImageData
#include <vtkImageData.h>
#include <vtkPolyDataToImageStencil.h>
#include <vtkImageStencil.h>
#include <vtkPointData.h>

// DisplaySSM Logic includes
#include "vtkSlicerDisplaySSMLogic.h"

#include <vtkNew.h>
#include "vtkMRMLModelNode.h"
#include "vtkMRMLModelDisplayNode.h"
#include "vtkMRMLScene.h"
#include "vtkMRMLLayoutNode.h"
#include "vtkMRMLDoubleArrayNode.h"
#include "vtkDoubleArray.h"
#include "vtkMRMLChartViewNode.h"
#include "vtkMRMLChartNode.h"

#include "vtkMRMLScalarVolumeNode.h"
#include "vtkMRMLScalarVolumeDisplayNode.h"
#include "vtkMRMLLabelMapVolumeDisplayNode.h"

#include "vtkMRMLSelectionNode.h"
#include "vtkMatrix4x4.h"
#include <vtkImageChangeInformation.h>

//
#include "qSlicerApplication.h"
#include "qSlicerLayoutManager.h"
#include "qMRMLThreeDView.h"
#include "qSlicerIO.h"
#include "qSlicerCoreIOManager.h"

#include <time.h>


typedef itk::Mesh<float, Dimensions> MeshType;
typedef vnl_vector<statismo::ScalarType> itkVectorType;

//-----------------------------------------------------------------------------
/// \ingroup Slicer_QtModules_ExtensionTemplate
class qSlicerDisplaySSMModuleWidgetPrivate: public Ui_qSlicerDisplaySSMModuleWidget
{
public:
  qSlicerDisplaySSMModuleWidgetPrivate();

  vtkPolyData* convertMeshToVtk(MeshType::Pointer meshToConvert, vtkPolyData *  m_PolyDataReturn);

  vtkImageData* convertPolyDataToImageData(vtkPolyData* inpuPolyData, double spacing[], double *bounds);

};

//-----------------------------------------------------------------------------
// qSlicerDisplaySSMModuleWidgetPrivate methods

//-----------------------------------------------------------------------------
qSlicerDisplaySSMModuleWidgetPrivate::qSlicerDisplaySSMModuleWidgetPrivate()
{
}

//-----------------------------------------------------------------------------
vtkPolyData* qSlicerDisplaySSMModuleWidgetPrivate::convertMeshToVtk(
  MeshType::Pointer meshToConvert, vtkPolyData *  m_PolyDataReturn)
{
  typedef MeshType::MeshTraits                      TriangleMeshTraits;
  typedef MeshType::PointType                       PointType;
  typedef MeshType::PointsContainer                 InputPointsContainer;
  typedef InputPointsContainer::Pointer             InputPointsContainerPointer;
  typedef InputPointsContainer::Iterator            InputPointsContainerIterator;
  typedef MeshType::CellType                        CellType;

  typedef MeshType::CellsContainerPointer           CellsContainerPointer;
  typedef MeshType::CellsContainerIterator          CellsContainerIterator;

  vtkSmartPointer< vtkPoints >   m_Points = vtkSmartPointer< vtkPoints >::New();
  vtkSmartPointer< vtkPolyData >   m_PolyData = vtkSmartPointer< vtkPolyData >::New();
  vtkCellArray * m_Polys = vtkCellArray::New();

  int numPoints =  meshToConvert->GetNumberOfPoints();

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
vtkImageData* qSlicerDisplaySSMModuleWidgetPrivate::convertPolyDataToImageData(
  vtkPolyData* inpuPolyData, double spacing[], double *bounds)
{
  vtkNew< vtkImageData >   whiteImage;
  inpuPolyData->GetBounds(bounds);
  whiteImage->SetSpacing(spacing);

  // compute dimensions
  int dim[3];
  for (int i = 0; i < 3; i++)
  {
    dim[i] = static_cast<int>(ceil((bounds[i * 2 + 1] - bounds[i * 2]) / spacing[i]));
  }

  whiteImage->SetDimensions(dim);
  whiteImage->SetExtent(0, dim[0] - 1, 0, dim[1] - 1, 0, dim[2] - 1);

  double origin[3];
  origin[0] = bounds[0] + spacing[0] / 2;
  origin[1] = bounds[2] + spacing[1] / 2;
  origin[2] = bounds[4] + spacing[2] / 2;
  whiteImage->SetOrigin(origin);
  whiteImage->SetScalarTypeToUnsignedChar();
  whiteImage->AllocateScalars();

  // fill the image with foreground voxels:
  unsigned char inval = 128; //255
  unsigned char outval = 0;
  vtkIdType count = whiteImage->GetNumberOfPoints();
  for (vtkIdType i = 0; i < count; ++i)
  {
    whiteImage->GetPointData()->GetScalars()->SetTuple1(i, inval);
  }

  // polygonal data --> image stencil:
  vtkNew< vtkPolyDataToImageStencil > pol2stenc;
  pol2stenc->SetInput(inpuPolyData);
  pol2stenc->SetOutputOrigin(origin);
  pol2stenc->SetOutputSpacing(spacing);
  pol2stenc->SetOutputWholeExtent(whiteImage->GetExtent());
  pol2stenc->Update();

  // cut the corresponding white image and set the background:
  vtkNew< vtkImageStencil > imgstenc;
  imgstenc->SetInput(whiteImage.GetPointer());
  imgstenc->SetStencil(pol2stenc->GetOutput());
  imgstenc->ReverseStencilOff();
  imgstenc->SetBackgroundValue(outval);
  imgstenc->Update();

  vtkImageData * outputImageData = imgstenc->GetOutput();
  outputImageData->Register(0);

  return outputImageData;
}



//-----------------------------------------------------------------------------
// qSlicerDisplaySSMModuleWidget methods

//-----------------------------------------------------------------------------
qSlicerDisplaySSMModuleWidget::qSlicerDisplaySSMModuleWidget(QWidget* _parent)
  : Superclass( _parent )
  , d_ptr( new qSlicerDisplaySSMModuleWidgetPrivate )
{
  itkModel = ItkStatisticalModelType::New();
  //vtkModel(new VtkStatisticalModelType);
}

//-----------------------------------------------------------------------------
qSlicerDisplaySSMModuleWidget::~qSlicerDisplaySSMModuleWidget()
{
  Q_D(qSlicerDisplaySSMModuleWidget);
  if (d->radioButtonVTK->isChecked()){
    vtkModel->Delete();
  }
}

//-----------------------------------------------------------------------------
void qSlicerDisplaySSMModuleWidget::setup()
{
  Q_D(qSlicerDisplaySSMModuleWidget);
  d->setupUi(this);
}

//-----------------------------------------------------------------------------
void qSlicerDisplaySSMModuleWidget::setMRMLScene(vtkMRMLScene* mrmlScene)
{
  this->qSlicerAbstractModuleWidget::setMRMLScene(mrmlScene);
}

//-----------------------------------------------------------------------------
void qSlicerDisplaySSMModuleWidget::onSelectInputModel()
{
  Q_D(qSlicerDisplaySSMModuleWidget);
  QString inputFile = QFileDialog::getOpenFileName(this, "Select input model", QString());
  d->modelNamePath->setText(inputFile);
}

//-----------------------------------------------------------------------------
void qSlicerDisplaySSMModuleWidget::applyModel()
{
  Q_D(qSlicerDisplaySSMModuleWidget);
  vtkPolyData* meanModel = vtkPolyData::New();
  unsigned int nbPrincipalComponent = 0;
  // Load the model
  std::string modelString = d->modelNamePath->text().toStdString();

  if (d->radioButtonVTK->isChecked()){
    try {
      vtkModel = VtkStatisticalModelType::Load(modelString.c_str());
      meanModel = vtkModel->DrawMean();
      nbPrincipalComponent = vtkModel->GetNumberOfPrincipalComponents();
    }
    catch (StatisticalModelException& e) {
      std::cout << "Exception occured while building the shape model" << std::endl;
      std::cout << e.what() << std::endl;
      return;
    }
  }
  if (d->radioButtonITK->isChecked()){
    try {
      itkModel->Load(modelString.c_str());
      nbPrincipalComponent = itkModel->GetNumberOfPrincipalComponents();
      //Calculate mean
      typedef ItkRepresenterType::MeshType TestType;
      TestType::Pointer meanDf = itkModel->DrawMean();
      meanModel = d->convertMeshToVtk(meanDf, meanModel);
      nameITK.push_back("Mean");
      indexNode=0;
    }
    catch (itk::ExceptionObject& o) {
      std::cout << "Exception occured while building the shape model" << std::endl;
      std::cout << o << std::endl;
      return ;
    }
  }

  if (!d->radioButtonITK->isChecked() && !d->radioButtonVTK->isChecked()){
      return;
    }

  // Display Model and Volume
  displayModelVolume(meanModel, "Mean");
  // Add mean model to the scene
  /*vtkSlicerDisplaySSMLogic* moduleLogic = vtkSlicerDisplaySSMLogic::New();
  moduleLogic->DisplaySampleModel(meanModel, this->mrmlScene());*/
  
  meanModel->Delete();

  // Display Eigen spectrum
  displayEigenSpectrum(nbPrincipalComponent);

  // Set the number of components for the pcSlider
  d->pcSlider->setMaximum(nbPrincipalComponent);
  d->pcSlider->setValue(1);
  d->stdSlider->setValue(0);

}

//-----------------------------------------------------------------------------
void qSlicerDisplaySSMModuleWidget::onSelect()
{
  /**** Display Model ****/
  using std::auto_ptr;
  Q_D(qSlicerDisplaySSMModuleWidget);
  vtkPolyData* samplePC = vtkPolyData::New();

  int pc = d->pcSlider->value()-1; // -1 because user can choose between 1 and max
  int std = d->stdSlider->value();

  //display Sample name with PC and Std values
  std::ostringstream ssSample;
  ssSample <<"samplePc"<<pc<<"Std"<<std;
  std::string sampleName = ssSample.str();

  //Get the model Display node
  vtkSmartPointer<vtkCollection> modelDisplayNodes = vtkSmartPointer<vtkCollection>::Take( this->mrmlScene()->GetNodesByClass("vtkMRMLModelDisplayNode") );
  vtkSmartPointer<vtkCollection> volumeNodes = vtkSmartPointer<vtkCollection>::Take( this->mrmlScene()->GetNodesByClass("vtkMRMLScalarVolumeNode") );
  
  if (d->radioButtonVTK->isChecked()){
    int nbPrincipalComponent = vtkModel->GetNumberOfPrincipalComponents();
    VectorType coefficients = VectorType::Zero(nbPrincipalComponent);
    coefficients(pc) = std;
    samplePC = vtkModel->DrawSample(coefficients);
    /*double prob = vtkModel->ComputeProbabilityOfDataset(samplePC);
    std::cout<<"prob = "<<prob<<std::endl;*/
  }
  if (d->radioButtonITK->isChecked()){

    typedef std::vector<std::string>::iterator ItemIterator;
    ItemIterator iter = std::find(nameITK.begin(), nameITK.end(), sampleName); //check if the name already exists
    size_t index = std::distance(nameITK.begin(), iter);

    //Set visibility off of the previous model
    vtkMRMLModelDisplayNode* modelViewNode = vtkMRMLModelDisplayNode::SafeDownCast( modelDisplayNodes->GetItemAsObject (3+indexNode) );
    modelViewNode->VisibilityOff();

    if(iter!=nameITK.end()){ //model already exists
      //Set Visibility on of the current model
      vtkMRMLModelDisplayNode* modelViewNode = vtkMRMLModelDisplayNode::SafeDownCast( modelDisplayNodes->GetItemAsObject (3+index) );
      modelViewNode->VisibilityOn();
      //Display the volume in the slice node
      vtkMRMLScalarVolumeNode* volumeNode = vtkMRMLScalarVolumeNode::SafeDownCast( volumeNodes->GetItemAsObject (index) );
  	  vtkSlicerApplicationLogic * appLogic = qSlicerCoreApplication::application()->applicationLogic();
  	  vtkMRMLSelectionNode * selectionNode = appLogic->GetSelectionNode();
  	  selectionNode->SetReferenceActiveVolumeID(volumeNode->GetID());
  	  appLogic->PropagateVolumeSelection();
  	  appLogic->FitSliceToAll();

    }else{
      nameITK.push_back(sampleName);

      int nbPrincipalComponent = itkModel->GetNumberOfPrincipalComponents();
      itkVectorType coefficients(nbPrincipalComponent,0.0); // initialize the vector to 0
      coefficients(pc) = std;
      //Calculate sample
      typedef ItkRepresenterType::MeshType TestType;
      TestType::Pointer itkSamplePC = itkModel->DrawSample(coefficients);

      /*double prob = itkModel->ComputeProbabilityOfDataset(itkSamplePC);
      std::cout<<"prob = "<<prob<<std::endl;*/

      samplePC=d->convertMeshToVtk(itkSamplePC, samplePC);
      displayModelVolume(samplePC, sampleName);
    }
    indexNode = index;
  }
  samplePC->Delete();

}

//-----------------------------------------------------------------------------
void qSlicerDisplaySSMModuleWidget::displayEigenSpectrum(unsigned int nbPrincipalComponent)
{
  Q_D(qSlicerDisplaySSMModuleWidget);
  std::vector<statismo::ScalarType> eigenvalueVector;
  double sumEigenvalue = 0;

  if (d->radioButtonVTK->isChecked()){
    VectorType eigenvalue = vtkModel->GetPCAVarianceVector();
    sumEigenvalue = eigenvalue.sum();
    for (unsigned int i=0;i<eigenvalue.rows();i++){
      eigenvalueVector.push_back(eigenvalue[i]);
    }
  }
  if (d->radioButtonITK->isChecked()){
    itkVectorType itkEigenvalue = itkModel->GetPCAVarianceVector();
    sumEigenvalue = itkEigenvalue.sum();
    for (unsigned int i=0;i<itkEigenvalue.size();i++){
      eigenvalueVector.push_back(itkEigenvalue(i));
    }
  }
  
  if (nbPrincipalComponent==eigenvalueVector.size()){
    std::cout<<"ok"<<std::endl;
  }else {
    std::cout<<"pas ok: nb Component = "<<nbPrincipalComponent<<"  size vector="<<eigenvalueVector.size()<<std::endl;
  }

  vtkSmartPointer<vtkCollection> layoutNodes = vtkSmartPointer<vtkCollection>::Take( this->mrmlScene()->GetNodesByClass("vtkMRMLLayoutNode") );
  layoutNodes->InitTraversal();
  vtkObject* layoutNodeVtkObject = layoutNodes->GetNextItemAsObject();
  vtkMRMLLayoutNode* layoutNode = vtkMRMLLayoutNode::SafeDownCast(layoutNodeVtkObject);
  if (!layoutNode)
  {
    std::cout<<"GetChartViewNode: Unable to get layout node!"<<std::endl;
  }
  layoutNode->SetViewArrangement( vtkMRMLLayoutNode::SlicerLayoutConventionalQuantitativeView );
  
  vtkSmartPointer<vtkCollection> chartViewNodes = vtkSmartPointer<vtkCollection>::Take( this->mrmlScene()->GetNodesByClass("vtkMRMLChartViewNode") );
  chartViewNodes->InitTraversal();
  vtkMRMLChartViewNode* chartViewNode = vtkMRMLChartViewNode::SafeDownCast( chartViewNodes->GetNextItemAsObject() );
  if (!chartViewNode)
  {
    std::cout<<"GetChartViewNode: Unable to get chart view node!"<<std::endl;
  }

  // Create an Array Node and add the eigen value
  vtkSmartPointer<vtkMRMLDoubleArrayNode> doubleArrayNode = vtkMRMLDoubleArrayNode::SafeDownCast(this->mrmlScene()->AddNode(vtkSmartPointer<vtkMRMLDoubleArrayNode>::New()));
  vtkDoubleArray* a = doubleArrayNode->GetArray();
  a->SetNumberOfTuples(nbPrincipalComponent);
  double cpt = 0;
  for (unsigned int i = 0;i<nbPrincipalComponent;i++){
    cpt = cpt + eigenvalueVector[i]*100/sumEigenvalue;
    a->SetComponent(i, 0, i);
    a->SetComponent(i, 1, cpt);
    a->SetComponent(i, 2, 0);
  }
  
  // Create a Chart Node.
  vtkNew<vtkMRMLChartNode> chartNode;
  chartNode->AddArray("EigenSpectrum", doubleArrayNode->GetID());
  this->mrmlScene()->AddNode(chartNode.GetPointer());
  
  // Set properties on the Chart
  chartNode->SetProperty("default", "title", "EigenSpectrum");
  chartNode->SetProperty("default", "xAxisLabel", "Principal components");
  chartNode->SetProperty("default", "yAxisLabel", "% of variation captured");
  chartNode->SetProperty("default", "type", "Line");
  chartNode->SetProperty("default", "showMarkers", "on");

  // Tell the Chart View which Chart to display
  chartViewNode->SetChartNodeID(chartNode->GetID());
}

//-----------------------------------------------------------------------------
void qSlicerDisplaySSMModuleWidget::displayModelVolume(vtkPolyData* modelToDisplay, std::string modelName)
{
  Q_D(qSlicerDisplaySSMModuleWidget);
  
  // Display model
  vtkNew<vtkMRMLModelDisplayNode> modelDisplayNode;
  this->mrmlScene()->AddNode(modelDisplayNode.GetPointer());

  vtkNew<vtkMRMLModelNode> meanNode;
  meanNode->SetAndObservePolyData(modelToDisplay);
  meanNode->SetAndObserveDisplayNodeID(modelDisplayNode->GetID());
  meanNode->SetName(modelName.c_str());
  this->mrmlScene()->AddNode(meanNode.GetPointer());
  
  double spacing[3]; // desired volume spacing
  spacing[0] = 0.24;
  spacing[1] = 0.24;
  spacing[2] = 0.6;

  double bounds[6];
  vtkSmartPointer<vtkImageData> meanImage;
  meanImage.TakeReference(d->convertPolyDataToImageData(modelToDisplay, spacing, &bounds[0]));

  //Create a displayable node and add to the scene
  vtkNew<vtkMRMLLabelMapVolumeDisplayNode> volumeDisplayNode;
  volumeDisplayNode->SetAndObserveColorNodeID("vtkMRMLColorTableNodeFileGenericAnatomyColors.txt");
  this->mrmlScene()->AddNode( volumeDisplayNode.GetPointer() );
  
  vtkNew<vtkImageChangeInformation> ici;
  ici->SetInput(meanImage);
  ici->SetOutputSpacing( 1, 1, 1 );
  ici->SetOutputOrigin( 0, 0, 0 );
  ici->Update();

  vtkNew<vtkMatrix4x4> matrix;
  matrix->Identity();
  matrix->SetElement(0,0,-1);
  matrix->SetElement(1,1,-1);
  
 //volume name with PC and Std values
  std::ostringstream ssVolume;
  ssVolume <<"volumePc"<<d->pcSlider->value()-1<<"Std"<<d->stdSlider->value();
  std::string volumeName = ssVolume.str();
  
  //Create a scalar volume node with the created volume
  vtkNew<vtkMRMLScalarVolumeNode> volumeNode;
  volumeNode->SetIJKToRASMatrix(matrix.GetPointer());
  volumeNode->SetAndObserveDisplayNodeID(volumeDisplayNode->GetID());
  volumeNode->SetAndObserveImageData(ici->GetOutput());
  volumeNode->SetOrigin(-bounds[0], -bounds[2], bounds[4]); //Compensate Slicer coordinates
  volumeNode->SetSpacing(spacing);
  volumeNode->SetLabelMap(true);
  volumeNode->SetName(volumeName.c_str());
  this->mrmlScene()->AddNode(volumeNode.GetPointer());

  // finally display the volume in the slice node
  vtkSlicerApplicationLogic * appLogic = qSlicerCoreApplication::application()->applicationLogic();
  vtkMRMLSelectionNode * selectionNode = appLogic->GetSelectionNode();
  //selectionNode->SetReferenceActiveLabelVolumeID(volumeNode->GetID());
  selectionNode->SetReferenceActiveVolumeID(volumeNode->GetID());
  appLogic->PropagateVolumeSelection();
  appLogic->FitSliceToAll();
}

