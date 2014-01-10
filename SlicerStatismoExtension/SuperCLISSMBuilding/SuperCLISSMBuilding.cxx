#include "SuperCLISSMBuildingCLP.h"

#include <vtksys/SystemTools.hxx>

#include "statismo/PCAModelBuilder.h"
#include "statismo/StatisticalModel.h"
#include "statismo/DataManager.h"

#include "Representers/VTK/vtkPolyDataRepresenter.h"

#include <vtkDirectory.h>
#include <vtkPolyData.h>
#include <vtkPolyDataWriter.h>

#include <iostream>
#include <memory>

//#include "vtkMRMLScene.h"
//#include "vtkMRMLModelNode.h"
//#include "vtkMRMLModelDisplayNode.h"

#include "vtkXMLPolyDataReader.h"
#include "vtkXMLPolyDataWriter.h"

using namespace statismo;
using std::auto_ptr;

typedef vtkPolyDataRepresenter RepresenterType;
typedef StatisticalModel<RepresenterType> StatisticalModelType;

vtkPolyData* loadVTKPolyData(const std::string& filename)
{
  vtkPolyDataReader* reader = vtkPolyDataReader::New();
  reader->SetFileName(filename.c_str());
  reader->Update();
  vtkPolyData* pd = vtkPolyData::New();
  pd->ShallowCopy(reader->GetOutput());
  return pd;
}

void displayModel(vtkPolyData* modelToDisplay, std::string modelDisplayName){
  // write the output
  std::string::size_type loc = modelDisplayName.find_last_of(".");
  std::string extension = modelDisplayName.substr(loc);
  if( loc == std::string::npos )
	{
	std::cout << "Failed to find an extension for " << modelDisplayName << std::endl;
	}
  if( extension == std::string(".vtk") )
    {
    vtkPolyDataWriter *pdWriter = vtkPolyDataWriter::New();
    pdWriter->SetFileName(modelDisplayName.c_str() );
    pdWriter->SetInput(modelToDisplay);
    pdWriter->Write();
    pdWriter->Update();
    pdWriter->Delete();
    }
  else if( extension == std::string(".vtp") )
    {
    vtkXMLPolyDataWriter *pdWriter = vtkXMLPolyDataWriter::New();
    pdWriter->SetIdTypeToInt32();
    pdWriter->SetFileName(modelDisplayName.c_str());
    pdWriter->SetInput(modelToDisplay);
    pdWriter->Write();
    pdWriter->Update();
    pdWriter->Delete();
    }
}


void modelBuilding(std::string refModel, std::string inputModelDir, std::string modelName, std::string meanModels){
	// All the statismo classes have to be parameterized with the RepresenterType.
	// For building a shape model with vtk, we use the vtkPolyDataRepresenter.
	typedef DataManager<RepresenterType> DataManagerType;
	typedef PCAModelBuilder<RepresenterType> ModelBuilderType;

	try {

	// We create a new representer object. For the vtkPolyDataRepresenter, we have to set a reference
	// and the alignmentType. The alignmenttype (which is here RIGID) determines how the dataset that we
	// will use will later be aligned to the reference.
    std::cout<<"ref model= "<< refModel<<std::endl;
	vtkPolyData* reference = loadVTKPolyData(inputModelDir + "/modelFitted0.vtk"); //why not a vtkpolydata? instead of vtp 
	auto_ptr<RepresenterType> representer(RepresenterType::Create(reference, RepresenterType::RIGID));

	// We create a datamanager and provide it with a pointer  to the representer
	auto_ptr<DataManagerType> dataManager(DataManagerType::Create(representer.get()));

    vtkDirectory* directory = vtkDirectory::New();
    int opened = directory->Open(inputModelDir.c_str());

    int numberOfFiles = directory->GetNumberOfFiles();
	// Now we add our data to the data manager
	// load the data and add it to the data manager. We take the first 17 hand shapes that we find in the data folder
	for (int i = 0; i < numberOfFiles; i++) {
      std::string fileString = inputModelDir;
      fileString += "/";
      fileString += directory->GetFile(i);
			std::cout << "file = " << fileString << std::endl;

      std::string ext = vtksys::SystemTools::GetFilenameLastExtension(fileString);
      std::cout << fileString << " extension: " << ext << std::endl;
      if (ext==".vtk"){
		  vtkPolyData* dataset = loadVTKPolyData(fileString);

		  // We provde the filename as a second argument.
		  // It will be written as metadata, and allows us to more easily figure out what we did later.
		  dataManager->AddDataset(dataset, fileString);

		  // it is save to delete the dataset after it was added, as the datamanager direclty copies it.
		  dataset->Delete();
      }
	}

	// To actually build a model, we need to create a model builder object.
	// Calling the build model with a list of samples from the data manager, returns a new model.
	// The second parameter to BuildNewModel is the variance of the noise on our data
	auto_ptr<ModelBuilderType> modelBuilder(ModelBuilderType::Create());

	auto_ptr<StatisticalModelType> model(modelBuilder->BuildNewModel(dataManager->GetSampleData(), 0.01));

	// Once we have built the model, we can save it to disk.
	model->Save(modelName);
	std::cout << "Successfully saved shape model as " << modelName << std::endl;

	reference->Delete();

  	// get the model mean
	vtkPolyData* mean = model->DrawMean();

    // write the output
	displayModel(mean, meanModels);
	}
	catch (StatisticalModelException& e) {
		std::cout << "Exception occured while building the shape model" << std::endl;
		std::cout << e.what() << std::endl;
	}

}


int main( int argc, char * argv[] )
{
  PARSE_ARGS;
  std::cout<<"skipModelBuilding= "<<skipModelBuilding<<std::endl;

  if(!skipModelBuilding){
    modelBuilding(referenceModel, inputModelDirectory, outputModelName, meanModel);
  }else{
    std::cout<<"display"<<std::endl;
    std::cout<<"outputModelName= "<<outputModelName<<std::endl;
    auto_ptr<StatisticalModelType> model(StatisticalModelType::Load(outputModelName));
	VectorType coefficients = VectorType::Zero(model->GetNumberOfPrincipalComponents());
    std::cout<<"model->GetNumberOfPrincipalComponents()= "<<model->GetNumberOfPrincipalComponents()<<std::endl;
	coefficients(PC) = standardVariation;
	vtkPolyData* samplePC = model->DrawSample(coefficients);
    displayModel(samplePC, modelSamplePC);
  }
  std::cout<<"OK"<<std::endl;
  return EXIT_SUCCESS;

}


