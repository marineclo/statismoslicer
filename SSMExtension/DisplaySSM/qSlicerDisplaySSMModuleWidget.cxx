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

//
#include "qSlicerApplication.h"
#include "qSlicerLayoutManager.h"
#include "qMRMLThreeDView.h"

#include <time.h>

const unsigned Dimensions = 3;
typedef itk::Mesh<float, Dimensions> MeshType;
typedef itk::MeshFileReader<MeshType> MeshFileReaderType;
typedef itk::MeshRepresenter<float, Dimensions> ItkRepresenterType;
typedef itk::StatisticalModel<ItkRepresenterType> ItkStatisticalModelType;
typedef vnl_vector<statismo::ScalarType> itkVectorType;

//-----------------------------------------------------------------------------
/// \ingroup Slicer_QtModules_ExtensionTemplate
class qSlicerDisplaySSMModuleWidgetPrivate: public Ui_qSlicerDisplaySSMModuleWidget
{
public:
  qSlicerDisplaySSMModuleWidgetPrivate();

  ItkStatisticalModelType* itkStatModel;
  void setItkStatModel(ItkStatisticalModelType* newItkStatModel);
  ItkStatisticalModelType* getItkStatModel();
  vtkPolyData* convertMeshToVtk(MeshType::Pointer meshToConvert, vtkPolyData *  m_PolyDataReturn);

};

//-----------------------------------------------------------------------------
// qSlicerDisplaySSMModuleWidgetPrivate methods

//-----------------------------------------------------------------------------
qSlicerDisplaySSMModuleWidgetPrivate::qSlicerDisplaySSMModuleWidgetPrivate()
{
}

void qSlicerDisplaySSMModuleWidgetPrivate::setItkStatModel(ItkStatisticalModelType* newItkStatModel)
{
    if (this->itkStatModel != 0)
    {
        this->itkStatModel->Delete();
        this->itkStatModel = 0;
        std::cout<<"itkStatModel not 0"<<std::endl;

    }
    this->itkStatModel = ItkStatisticalModelType::New();
    this->itkStatModel = newItkStatModel;
}

ItkStatisticalModelType* qSlicerDisplaySSMModuleWidgetPrivate::getItkStatModel()
{
  ItkStatisticalModelType* newItkStatModel = ItkStatisticalModelType::New();
  return newItkStatModel = this->itkStatModel;
}

vtkPolyData* qSlicerDisplaySSMModuleWidgetPrivate::convertMeshToVtk(MeshType::Pointer meshToConvert, vtkPolyData *  m_PolyDataReturn)
{
  typedef typename MeshType::MeshTraits                      TriangleMeshTraits;
  typedef typename MeshType::PointType                       PointType;
  typedef typename MeshType::PointsContainer                 InputPointsContainer;
  typedef typename InputPointsContainer::Pointer             InputPointsContainerPointer;
  typedef typename InputPointsContainer::Iterator            InputPointsContainerIterator;
  typedef typename MeshType::CellType                        CellType;

  typedef typename MeshType::CellsContainerPointer           CellsContainerPointer;
  typedef typename MeshType::CellsContainerIterator          CellsContainerIterator;

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
    typename CellType::PointIdIterator pointIt = nextCell->PointIdsBegin();
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
// qSlicerDisplaySSMModuleWidget methods

//-----------------------------------------------------------------------------
qSlicerDisplaySSMModuleWidget::qSlicerDisplaySSMModuleWidget(QWidget* _parent)
  : Superclass( _parent )
  , d_ptr( new qSlicerDisplaySSMModuleWidgetPrivate )
{
}

//-----------------------------------------------------------------------------
qSlicerDisplaySSMModuleWidget::~qSlicerDisplaySSMModuleWidget()
{
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
  // Load the model
  std::string modelString = inputFile.toStdString();
  try {
    ItkStatisticalModelType* modelITK = ItkStatisticalModelType::New();
    modelITK->Load(modelString.c_str());
    //d->setItkStatModel(modelITK);
    //d->setItkStatModel(d->itkStatModel->Load(modelString.c_str()));

    // Display Eigen spectrum
    itkVectorType itkEigenvalue = modelITK->GetPCAVarianceVector(); 
    unsigned int nbPrincipalComponent = modelITK->GetNumberOfPrincipalComponents();

    // Set the number of components for the pcSlider
    //d->pcSlider->setMaximum(d->itkStatModel->GetNumberOfPrincipalComponents());
    d->pcSlider->setMaximum(nbPrincipalComponent);
     std::cout<<"test5"<<std::endl;
    
    //Calculate mean
    typedef typename ItkRepresenterType::MeshType TestType;
    typename TestType::Pointer meanDf = modelITK->DrawMean();

    vtkPolyData* meanModel = vtkPolyData::New();
    meanModel=d->convertMeshToVtk(meanDf, meanModel);

    std::cout<<"test6"<<std::endl;
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

    displayEigenSpectrum(itkEigenvalue, nbPrincipalComponent);

  }
 catch (itk::ExceptionObject& o) {
  std::cout << "Exception occured while building the shape model" << std::endl;
  std::cout << o << std::endl;
 }
}

//-----------------------------------------------------------------------------
void qSlicerDisplaySSMModuleWidget::onSelect()
{
  /**** ITK Model ****/
  using std::auto_ptr;
  Q_D(qSlicerDisplaySSMModuleWidget);
  vtkPolyData* samplePC = vtkPolyData::New();
  // Get the model name selected by the user
  //d->itkStatModel = d->getItkStatModel();

  QString inputFile = d->modelNamePath->text();
  // Load the model
  std::string modelString = inputFile.toStdString();
  ItkStatisticalModelType* modelITKw = ItkStatisticalModelType::New();
  modelITKw->Load(modelString.c_str());

  int nbPrincipalComponent = modelITKw->GetNumberOfPrincipalComponents();
  itkVectorType coefficients(nbPrincipalComponent,0.0); // set the vector to 0
  int pc = d->pcSlider->value()-1; // -1 because user can choose between 1 and max
  int std = d->stdSlider->value();
  coefficients(pc) = std;
  //MeshType::Pointer itkSamplePC = d->itkStatModel->DrawSample(coefficients);
  //MeshType::Pointer mesh = d->itkStatModel->DrawSample(coefficients);
 
  //Calculate sample
  typedef typename ItkRepresenterType::MeshType TestType;
  typename TestType::Pointer itkSamplePC = modelITKw->DrawSample(coefficients);

  //double prob = d->itkStatModel->ComputeProbabilityOfDataset(itkSamplePC);
  /*double prob = modelITKw->ComputeProbabilityOfDataset(itkSamplePC);
  std::cout<<"prob = "<<prob<<std::endl;*/

  samplePC=d->convertMeshToVtk(itkSamplePC, samplePC);

  // Add polydata sample to the scene
  vtkNew<vtkMRMLModelDisplayNode> sampleDisplayNode;
  srand(time(NULL)); // initialisation de rand
  
  double r = (rand()/(double)RAND_MAX );
  double g = (rand()/(double)RAND_MAX );
  double b = (rand()/(double)RAND_MAX );
  std::cout<<"r = "<<r<<"g= "<<g<<" b= "<<b<<std::endl;
  sampleDisplayNode->SetColor(r,g,b); // 
  //sampleDisplayNode->SetColor(0.5,1,0); // green
  this->mrmlScene()->AddNode(sampleDisplayNode.GetPointer());

  vtkNew<vtkMRMLModelNode> sampleNode;
  sampleNode->SetAndObservePolyData(samplePC);
  sampleNode->SetAndObserveDisplayNodeID(sampleDisplayNode->GetID());

  std::ostringstream ssSample;
		ssSample <<"SamplePc"<<pc<<"Std"<<std;
		std::string sampleName = ssSample.str();

  sampleNode->SetName(sampleName.c_str());
  this->mrmlScene()->AddNode(sampleNode.GetPointer());

  samplePC->Delete();

  /*vtkSlicerDisplaySSMLogic* moduleLogic = vtkSlicerDisplaySSMLogic::New();
  moduleLogic->DisplaySampleModel(samplePC, this->mrmlScene());*/
  
}

//-----------------------------------------------------------------------------
void qSlicerDisplaySSMModuleWidget::displayEigenSpectrum(itkVectorType itkEigenvalue, unsigned int nbPrincipalComponent)
{
  std::vector<statismo::ScalarType> eigenvalueVector;

  double sumEigenvalue = 0;
  sumEigenvalue = itkEigenvalue.sum();
  for (unsigned int i=0;i<itkEigenvalue.size();i++){
    eigenvalueVector.push_back(itkEigenvalue(i));
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
  
  std::cout<<"test"<<std::endl;
  // Create a Chart Node.
  vtkNew<vtkMRMLChartNode> chartNode;
  chartNode->AddArray("EigenSpectrum", doubleArrayNode->GetID());
  this->mrmlScene()->AddNode(chartNode.GetPointer());

  std::cout<<"test2"<<std::endl;
  
  // Set properties on the Chart
  chartNode->SetProperty("default", "title", "EigenSpectrum");
  chartNode->SetProperty("default", "xAxisLabel", "Principal components");
  chartNode->SetProperty("default", "yAxisLabel", "Eigen value");
  chartNode->SetProperty("default", "type", "Line");
  chartNode->SetProperty("default", "showMarkers", "on");
  std::cout<<"test3"<<std::endl;
  // Tell the Chart View which Chart to display
  chartViewNode->SetChartNodeID(chartNode->GetID());

  std::cout<<"test4"<<std::endl;

}

