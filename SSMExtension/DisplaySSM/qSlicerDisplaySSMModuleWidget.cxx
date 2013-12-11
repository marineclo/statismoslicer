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
#include "itkMeshFileWriter.h"
#include "itkMeshFileReader.h"

#include <vtkPolyDataReader.h>

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

//Compare PolyData
#include "vtkPointLocator.h"
#include "vtkMath.h"

// DisplaySSM Logic includes
#include "vtkSlicerDisplaySSMLogic.h"

#include <vtkNew.h>
#include "vtkMRMLModelNode.h"
#include "vtkMRMLModelDisplayNode.h"
#include "vtkMRMLScene.h"
#include "vtkMRMLLayoutNode.h"
#include "vtkMRMLDoubleArrayNode.h"
#include <vtkDoubleArray.h>
#include "vtkMRMLChartViewNode.h"
#include "vtkMRMLChartNode.h"

#include "vtkMRMLScalarVolumeNode.h"
#include "vtkMRMLScalarVolumeDisplayNode.h"
#include "vtkMRMLLabelMapVolumeDisplayNode.h"

#include "vtkMRMLSelectionNode.h"
#include "vtkMatrix4x4.h"
#include "vtkImageReslice.h"

#include "vtkMRMLColorNode.h"
#include "vtkMRMLColorTableNode.h"
#include "vtkMRMLColorTableStorageNode.h"
#include "vtkLookupTable.h"

#include <vtkScalarBarActor.h>
#include <vtkScalarBarWidget.h>


//
#include "qSlicerApplication.h"
#include "qSlicerLayoutManager.h"
#include "qMRMLThreeDView.h"
#include "qSlicerIO.h"
#include "qSlicerCoreIOManager.h"

#include <vtkMetaImageWriter.h>

#include <time.h>

//const unsigned Dimensions = 3;
typedef itk::Mesh<float, Dimensions> MeshType;
typedef itk::MeshFileReader<MeshType> MeshFileReaderType;
//typedef itk::MeshRepresenter<float, Dimensions> ItkRepresenterType;
//typedef itk::StatisticalModel<ItkRepresenterType> ItkStatisticalModelType;
typedef vnl_vector<statismo::ScalarType> itkVectorType;

//-----------------------------------------------------------------------------
/// \ingroup Slicer_QtModules_ExtensionTemplate
class qSlicerDisplaySSMModuleWidgetPrivate: public Ui_qSlicerDisplaySSMModuleWidget
{
public:
  qSlicerDisplaySSMModuleWidgetPrivate();

  vtkPolyData* convertMeshToVtk(MeshType::Pointer meshToConvert, vtkPolyData *  m_PolyDataReturn);

  vtkImageData* convertPolyDataToImageData(vtkSmartPointer<vtkPolyData> inpuPolyData, vtkImageData* meanImageReturn,  double *origin, double spacing[], double *bounds);

};

//-----------------------------------------------------------------------------
// qSlicerDisplaySSMModuleWidgetPrivate methods

//-----------------------------------------------------------------------------
qSlicerDisplaySSMModuleWidgetPrivate::qSlicerDisplaySSMModuleWidgetPrivate()
{
}

vtkPolyData* qSlicerDisplaySSMModuleWidgetPrivate::convertMeshToVtk(MeshType::Pointer meshToConvert, vtkPolyData *  m_PolyDataReturn)
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

vtkImageData* qSlicerDisplaySSMModuleWidgetPrivate::convertPolyDataToImageData(vtkSmartPointer< vtkPolyData > inpuPolyData, vtkImageData* meanImageReturn, double *origin, double spacing[], double *bounds){

  vtkSmartPointer< vtkImageData >   whiteImage = vtkSmartPointer< vtkImageData >::New();
  //double bounds[6];
  inpuPolyData->GetBounds(bounds);

  whiteImage->SetSpacing(spacing);

  // compute dimensions
  int dim[3];
  for (int i = 0; i < 3; i++)
    {
    std::cout<<"bounds"<<i<<"="<<bounds[i]<<"bounds"<<i+3<<"="<<bounds[i+3]<<std::endl;
    dim[i] = static_cast<int>(ceil((bounds[i * 2 + 1] - bounds[i * 2]) / spacing[i]));
    }

  whiteImage->SetDimensions(dim);
  whiteImage->SetExtent(0, dim[0] - 1, 0, dim[1] - 1, 0, dim[2] - 1);

  //double origin[3];
  origin[0] = bounds[0] + spacing[0] / 2;
  origin[1] = bounds[2] + spacing[1] / 2;
  origin[2] = bounds[4] + spacing[2] / 2;
  std::cout<<"origin[0]= " <<origin[0]<<std::endl;
  std::cout<<"origin[1]= " <<origin[1]<<std::endl;
  std::cout<<"origin[2]= " <<origin[2]<<std::endl;
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
  vtkSmartPointer< vtkPolyDataToImageStencil > pol2stenc = vtkSmartPointer< vtkPolyDataToImageStencil >::New();

  pol2stenc->SetInput(inpuPolyData);

  pol2stenc->SetOutputOrigin(origin);
  pol2stenc->SetOutputSpacing(spacing);
  pol2stenc->SetOutputWholeExtent(whiteImage->GetExtent());
  pol2stenc->Update();

  // cut the corresponding white image and set the background:
  vtkSmartPointer< vtkImageStencil > imgstenc = vtkSmartPointer< vtkImageStencil >::New();

  imgstenc->SetInput(whiteImage);
  imgstenc->SetStencil(pol2stenc->GetOutput());

  imgstenc->ReverseStencilOff();
  imgstenc->SetBackgroundValue(outval);
  imgstenc->Update();

  meanImageReturn->ShallowCopy(imgstenc->GetOutput());

  return meanImageReturn;
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
  //vtkPolyData* meanModel = vtkPolyData::New();
  //unsigned int nbPrincipalComponent = 0;
  QString inputFile = QFileDialog::getOpenFileName(this, "Select input model", QString());
  d->modelNamePath->setText(inputFile);
}

void qSlicerDisplaySSMModuleWidget::applyModel()
{
  Q_D(qSlicerDisplaySSMModuleWidget);
  vtkPolyData* meanModel = vtkPolyData::New();
  unsigned int nbPrincipalComponent = 0;
  // Load the model
  //std::string modelString = inputFile.toStdString();
  std::string modelString = d->modelNamePath->text().toStdString();

  if (d->radioButtonVTK->isChecked()){
    try {
      //vtkModel(new VtkStatisticalModelType::Load(modelString.c_str()));
      //vtkModel->Load(modelString.c_str());
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

      meanModel=d->convertMeshToVtk(meanDf, meanModel);

      std::cout<<"test6"<<std::endl;

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


  // Add mean model to the scene
  /*vtkSlicerDisplaySSMLogic* moduleLogic = vtkSlicerDisplaySSMLogic::New();
  moduleLogic->DisplaySampleModel(meanModel, this->mrmlScene());*/

  std::cout<<"test7"<<std::endl;

  vtkNew<vtkMRMLModelDisplayNode> modelDisplayNode;
  this->mrmlScene()->AddNode(modelDisplayNode.GetPointer());

  vtkNew<vtkMRMLModelNode> meanNode;
  meanNode->SetAndObservePolyData(meanModel);
  meanNode->SetAndObserveDisplayNodeID(modelDisplayNode->GetID());
  meanNode->SetName("Mean");
  this->mrmlScene()->AddNode(meanNode.GetPointer());

  double spacing[3]; // desired volume spacing
  spacing[0] = 0.24;
  spacing[1] = 0.24;
  spacing[2] = 0.6;


  double origin[3];
  double bounds[6];
  vtkImageData* meanImage = vtkImageData::New();
  meanImage = d->convertPolyDataToImageData(meanModel, meanImage, &origin[0], spacing, &bounds[0]);
  std::cout<<"test5="<<origin[0]<<" "<<origin[1]<<" "<<origin[2]<<std::endl;

  //int* extentsMean = meanImage->GetExtent();
  //std::cout<<"extentsMean "<<extentsMean[0]<<" "<<extentsMean[1]<<" "<<extentsMean[2]<<" "<<extentsMean[3]<<" "<<extentsMean[4]<<" "<<extentsMean[5]<<std::endl;

  /*vtkMetaImageWriter* writer = vtkMetaImageWriter::New();
  //writer->SetFileName("/home/marine/Documents/SSM/Tools/PolyDataToImageData/meanVolume.mhd");
  writer->SetFileName("/Users/Marine/EPFL/Documents/Tools/ssmtool/PolyDataToImageData/meanVolume.mhd");
  writer->SetInput(meanImage);
  writer->Write();*/


  /*qSlicerCoreIOManager* coreIOManager = qSlicerCoreApplication::application()->coreIOManager();
  qSlicerIO::IOProperties fileParameters;
  fileParameters["filename"] = "/home/marine/Documents/SSM/Tools/PolyDataToImageData/meanVolume.mhd";
  vtkMRMLNode* volumeNode = coreIOManager->loadNodesAndGetFirst("VolumeFile", fileParameters);*/


  //Create a displayable node and add to the scene
  vtkNew<vtkMRMLLabelMapVolumeDisplayNode> volumeDisplayNode;
  volumeDisplayNode->SetAndObserveColorNodeID("vtkMRMLColorTableNodeFileGenericAnatomyColors.txt");
  this->mrmlScene()->AddNode( volumeDisplayNode.GetPointer() );
  //meanImage->SetOrigin(-bounds[0],-bounds[2],bounds[4]);
  meanImage->SetOrigin(0, 0, 0);
  double* origin2 = meanImage->GetOrigin();
  std::cout<<"test7="<<origin2[0]<<" "<<origin2[1]<<" "<<origin2[2]<<std::endl;

  //Create a scalar volume node with the created volume
  vtkNew<vtkMRMLScalarVolumeNode> volumeNode;

  vtkMatrix4x4* matrix = vtkMatrix4x4::New();
  matrix->Identity();
  matrix->SetElement(0,0,-1);
  matrix->SetElement(1,1,-1);

  volumeNode->SetIJKToRASMatrix(matrix);
  volumeNode->SetAndObserveDisplayNodeID(volumeDisplayNode->GetID());
  volumeNode->SetAndObserveImageData(meanImage);
  volumeNode->SetOrigin(-bounds[0], -bounds[2], bounds[4]); //Compensate Slicer coordinates
  //volumeNode->SetOrigin(14,36,bounds[4]); //Compensate Slicer coordinates
  volumeNode->SetSpacing(spacing);
  volumeNode->SetLabelMap(true);
  this->mrmlScene()->AddNode(volumeNode.GetPointer());

  // finally display the volume in the slice node
  vtkSlicerApplicationLogic * appLogic = qSlicerCoreApplication::application()->applicationLogic();

  vtkMRMLSelectionNode * selectionNode = appLogic->GetSelectionNode();
  //selectionNode->SetReferenceActiveLabelVolumeID(volumeNode->GetID());
  selectionNode->SetReferenceActiveVolumeID(volumeNode->GetID());
  appLogic->PropagateVolumeSelection();
  appLogic->FitSliceToAll();

  /*qSlicerApplication * app = qSlicerApplication::application();
  if (!app)
  {
    std::cout<<"error application"<<std::endl;
  }
  qSlicerLayoutManager * layoutManager = app->layoutManager();
  if (!layoutManager)
  {
    std::cout<<"error layout"<<std::endl;
  }
  qMRMLThreeDView* threeDView = layoutManager->threeDWidget(0)->threeDView();
  threeDView->resetFocalPoint();*/

  meanModel->Delete();
  matrix->Delete();
  //writer->Delete();
  meanImage->Delete();

   // Display Eigen spectrum
  displayEigenSpectrum(nbPrincipalComponent);

  // Set the number of components for the pcSlider
  d->pcSlider->setMaximum(nbPrincipalComponent);
  d->pcSlider->setValue(1);
  d->stdSlider->setValue(0);
  std::cout<<"test5"<<std::endl;

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
  ssSample <<"SamplePc"<<pc<<"Std"<<std;
  std::string sampleName = ssSample.str();

  bool alreadyCompute = false; //boolean to check if the model aready exist

  //Get the model Display node
  vtkSmartPointer<vtkCollection> modelDisplayNodes = vtkSmartPointer<vtkCollection>::Take( this->mrmlScene()->GetNodesByClass("vtkMRMLModelDisplayNode") );
  std::cout<<"Nb Collection= "<<modelDisplayNodes->GetNumberOfItems ()<<std::endl;

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

    if (nameITK.size()!=0){
      vtkMRMLModelDisplayNode* modelViewNode = vtkMRMLModelDisplayNode::SafeDownCast( modelDisplayNodes->GetItemAsObject (4+indexNode) );
      modelViewNode->VisibilityOff();
    }
    else{
      vtkMRMLModelDisplayNode* modelViewNode = vtkMRMLModelDisplayNode::SafeDownCast( modelDisplayNodes->GetItemAsObject (modelDisplayNodes->GetNumberOfItems ()-1) );
      modelViewNode->VisibilityOff();
    }

    if(iter!=nameITK.end()){ //model already exists
      alreadyCompute = true;
      //Set Visibility on of the current model
      vtkMRMLModelDisplayNode* modelViewNode = vtkMRMLModelDisplayNode::SafeDownCast( modelDisplayNodes->GetItemAsObject (4+index) );
      modelViewNode->VisibilityOn();

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
    }
    indexNode = index;
  }
  
  if (!alreadyCompute){

    // Add polydata sample to the scene
    vtkNew<vtkMRMLModelDisplayNode> sampleDisplayNode;
    srand(time(NULL)); // initialisation de rand

    double r = (rand()/(double)RAND_MAX ); //random color between 0 and 1
    double g = (rand()/(double)RAND_MAX );
    double b = (rand()/(double)RAND_MAX );
    sampleDisplayNode->SetColor(r,g,b); // random color
    this->mrmlScene()->AddNode(sampleDisplayNode.GetPointer());

    vtkNew<vtkMRMLModelNode> sampleNode;
    sampleNode->SetAndObservePolyData(samplePC);
    sampleNode->SetAndObserveDisplayNodeID(sampleDisplayNode->GetID());

    sampleNode->SetName(sampleName.c_str());
    this->mrmlScene()->AddNode(sampleNode.GetPointer());

    /*vtkSlicerDisplaySSMLogic* moduleLogic = vtkSlicerDisplaySSMLogic::New();
    moduleLogic->DisplaySampleModel(samplePC, this->mrmlScene());*/
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
  chartNode->SetProperty("default", "yAxisLabel", "Eigen value");
  chartNode->SetProperty("default", "type", "Line");
  chartNode->SetProperty("default", "showMarkers", "on");

  // Tell the Chart View which Chart to display
  chartViewNode->SetChartNodeID(chartNode->GetID());

}

//-----------------------------------------------------------------------------
void qSlicerDisplaySSMModuleWidget::selectPolyData1()
{
  Q_D(qSlicerDisplaySSMModuleWidget);
  QString inputFile = QFileDialog::getOpenFileName(this, "Select input polyData1", QString());
  d->polyData1Name->setText(inputFile);
}

//-----------------------------------------------------------------------------
void qSlicerDisplaySSMModuleWidget::selectPolyData2()
{
  Q_D(qSlicerDisplaySSMModuleWidget);
  QString inputFile = QFileDialog::getOpenFileName(this, "Select input model", QString());
  d->polyData2Name->setText(inputFile);
}

//-----------------------------------------------------------------------------
void qSlicerDisplaySSMModuleWidget::comparePolyData()
{
  Q_D(qSlicerDisplaySSMModuleWidget);
  std::string polyData1String = d->polyData1Name->text().toStdString();
  std::string polyData2String = d->polyData2Name->text().toStdString();

  vtkNew<vtkPolyDataReader> reader1;
  reader1->SetFileName(polyData1String.c_str());
  reader1->Update();
  vtkSmartPointer< vtkPolyData > polyData1 = reader1->GetOutput();

  vtkNew<vtkPolyDataReader> reader2;
  reader2->SetFileName(polyData2String.c_str());
  reader2->Update();
  vtkSmartPointer< vtkPolyData >  polyData2 = reader2->GetOutput();

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
  polyData1->GetPointData()->SetScalars(scalars);
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
  sampleNode->SetAndObservePolyData(polyData1);
  sampleNode->SetAndObserveDisplayNodeID(sampleDisplayNode->GetID());

  sampleNode->SetName(d->outputPolyDataName->text().toStdString().c_str());
  this->mrmlScene()->AddNode(sampleNode.GetPointer());

 /* vtkNew<vtkMRMLColorNode> colorNode;
  vtkScalarBarWidget* ScalarBarWidget = vtkScalarBarWidget::New();
  vtkScalarBarActor* ScalarBarActor = vtkScalarBarActor::New();
  ScalarBarWidget->SetScalarBarActor(ScalarBarActor);
  ScalarBarWidget->GetScalarBarActor()->SetOrientationToVertical();
  ScalarBarWidget->GetScalarBarActor()->SetNumberOfLabels(4);
  ScalarBarWidget->GetScalarBarActor()->SetMaximumNumberOfColors(256);
  ScalarBarWidget->GetScalarBarActor()->SetTitle("Error");
  ScalarBarWidget->GetScalarBarActor()->SetLabelFormat(" %s");

  // it's a 2d actor, position it in screen space by percentages
  ScalarBarWidget->GetScalarBarActor()->SetPosition(0.1, 0.1);
  ScalarBarWidget->GetScalarBarActor()->SetWidth(0.1);
  ScalarBarWidget->GetScalarBarActor()->SetHeight(0.8);
  ScalarBarWidget->GetScalarBarActor()->SetLookupTable();

  qSlicerApplication * app = qSlicerApplication::application();
  if (app && app->layoutManager())
  {
    qMRMLThreeDView* threeDView = app->layoutManager()->threeDWidget(0)->threeDView();
    vtkRenderer* activeRenderer = app->layoutManager()->activeThreeDRenderer();
    if (activeRenderer)
    {
      ScalarBarWidget->SetInteractor(activeRenderer->GetRenderWindow()->GetInteractor());
    }

  }
  //ScalarBarWidget->SetEnabled(true);*/


}
