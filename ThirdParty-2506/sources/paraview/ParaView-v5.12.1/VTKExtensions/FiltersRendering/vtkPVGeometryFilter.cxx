// SPDX-FileCopyrightText: Copyright (c) Kitware Inc.
// SPDX-License-Identifier: BSD-3-Clause

#include "vtkPVGeometryFilter.h"

#include "vtkAMRInformation.h"
#include "vtkAlgorithmOutput.h"
#include "vtkAppendPolyData.h"
#include "vtkCallbackCommand.h"
#include "vtkCellArray.h"
#include "vtkCellArrayIterator.h"
#include "vtkCellData.h"
#include "vtkCellGrid.h"
#include "vtkCellTypes.h"
#include "vtkCommand.h"
#include "vtkCompositeDataPipeline.h"
#include "vtkCompositeDataSet.h"
#include "vtkDataAssembly.h"
#include "vtkDataObjectTreeIterator.h"
#include "vtkExplicitStructuredGrid.h"
#include "vtkExplicitStructuredGridSurfaceFilter.h"
#include "vtkFeatureEdges.h"
#include "vtkFloatArray.h"
#include "vtkGarbageCollector.h"
#include "vtkGenericDataSet.h"
#include "vtkGenericGeometryFilter.h"
#include "vtkGeometryFilter.h"
#include "vtkHyperTreeGrid.h"
#include "vtkHyperTreeGridGeometry.h"
#include "vtkImageData.h"
#include "vtkInformation.h"
#include "vtkInformationIntegerVectorKey.h"
#include "vtkInformationVector.h"
#include "vtkMath.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkMultiPieceDataSet.h"
#include "vtkMultiProcessController.h"
#include "vtkNew.h"
#include "vtkObjectFactory.h"
#include "vtkOutlineSource.h"
#include "vtkOverlappingAMR.h"
#include "vtkPVTrivialProducer.h"
#include "vtkPartitionedDataSetCollection.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkPolygon.h"
#include "vtkRecoverGeometryWireframe.h"
#include "vtkRectilinearGrid.h"
#include "vtkRectilinearGridOutlineFilter.h"
#include "vtkSmartPointer.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkStringArray.h"
#include "vtkStructuredGrid.h"
#include "vtkStructuredGridOutlineFilter.h"
#include "vtkTimerLog.h"
#include "vtkTriangleFilter.h"
#include "vtkUniformGrid.h"
#include "vtkUnsignedIntArray.h"
#include "vtkUnstructuredGrid.h"
#include "vtkUnstructuredGridGeometryFilter.h"

#include <cassert>
#include <cmath>
#include <string>
#include <vector>

namespace details
{
static constexpr const char* ORIGINAL_FACE_IDS = "RecoverWireframeOriginalFaceIds";
}

template <typename T>
void GetValidWholeExtent(T* ds, const int wholeExt[6], int validWholeExt[6])
{
  if (wholeExt[0] <= wholeExt[1] && wholeExt[2] <= wholeExt[3] && wholeExt[4] <= wholeExt[5])
  {
    std::copy(wholeExt, wholeExt + 6, validWholeExt);
  }
  else
  {
    ds->GetExtent(validWholeExt);
  }
}

vtkStandardNewMacro(vtkPVGeometryFilter);
vtkCxxSetObjectMacro(vtkPVGeometryFilter, Controller, vtkMultiProcessController);
vtkInformationKeyMacro(vtkPVGeometryFilter, POINT_OFFSETS, IntegerVector);
vtkInformationKeyMacro(vtkPVGeometryFilter, VERTS_OFFSETS, IntegerVector);
vtkInformationKeyMacro(vtkPVGeometryFilter, LINES_OFFSETS, IntegerVector);
vtkInformationKeyMacro(vtkPVGeometryFilter, POLYS_OFFSETS, IntegerVector);
vtkInformationKeyMacro(vtkPVGeometryFilter, STRIPS_OFFSETS, IntegerVector);
class vtkPVGeometryFilter::BoundsReductionOperation : public vtkCommunicator::Operation
{
public:
  // Subclasses must overload this method, which performs the actual
  // operations.  The methods should first do a reinterpret_cast of the arrays
  // to the type suggested by \c datatype (which will be one of the VTK type
  // identifiers like VTK_INT, etc.).  Both arrays are considered top be
  // length entries.  The method should perform the operation A*B (where * is
  // a placeholder for whatever operation is actually performed) and store the
  // result in B.  The operation is assumed to be associative.  Commutativity
  // is specified by the Commutative method.
  void Function(const void* A, void* B, vtkIdType length, int datatype) override
  {
    assert((datatype == VTK_DOUBLE) && (length == 6));
    (void)datatype;
    (void)length;
    const double* bdsA = reinterpret_cast<const double*>(A);
    double* bdsB = reinterpret_cast<double*>(B);
    if (bdsA[0] < bdsB[0])
    {
      bdsB[0] = bdsA[0];
    }
    if (bdsA[1] > bdsB[1])
    {
      bdsB[1] = bdsA[1];
    }
    if (bdsA[2] < bdsB[2])
    {
      bdsB[2] = bdsA[2];
    }
    if (bdsA[3] > bdsB[3])
    {
      bdsB[3] = bdsA[3];
    }
    if (bdsA[4] < bdsB[4])
    {
      bdsB[4] = bdsA[4];
    }
    if (bdsA[5] > bdsB[5])
    {
      bdsB[5] = bdsA[5];
    }
  }

  // Description:
  // Subclasses override this method to specify whether their operation
  // is commutative.  It should return 1 if commutative or 0 if not.
  int Commutative() override { return 1; }
};

//----------------------------------------------------------------------------
vtkPVGeometryFilter::vtkPVGeometryFilter()
{
  this->OutlineFlag = 0;
  this->UseOutline = 1;
  this->GenerateFeatureEdges = false;
  this->BlockColorsDistinctValues = 7;
  // generating cell normals by default really slows down paraview
  // it is especially noticeable with the OpenGL2 backend.  Leaving
  // it on for the old backend as some tests rely on the cell normals
  // to be there as they use them for other purposes/etc.
  this->GenerateCellNormals = 0;
  this->Triangulate = false;
  this->NonlinearSubdivisionLevel = 1;

  this->GeometryFilter = vtkGeometryFilter::New();
  // we're prepping geometry for rendering
  // fast mode might generate wrong results because of certain assumptions that don't always hold
  // for that reason, the default is to have fast mode off
  this->GeometryFilter->SetFastMode(false);
  this->GenericGeometryFilter = vtkGenericGeometryFilter::New();
  this->UnstructuredGridGeometryFilter = vtkUnstructuredGridGeometryFilter::New();
  this->RecoverWireframeFilter = vtkRecoverGeometryWireframe::New();
  this->FeatureEdgesFilter = vtkFeatureEdges::New();

  // Setup a callback for the internal readers to report progress.
  this->GeometryFilter->AddObserver(
    vtkCommand::ProgressEvent, this, &vtkPVGeometryFilter::HandleGeometryFilterProgress);
  this->GenericGeometryFilter->AddObserver(
    vtkCommand::ProgressEvent, this, &vtkPVGeometryFilter::HandleGeometryFilterProgress);
  this->UnstructuredGridGeometryFilter->AddObserver(
    vtkCommand::ProgressEvent, this, &vtkPVGeometryFilter::HandleGeometryFilterProgress);
  this->RecoverWireframeFilter->AddObserver(
    vtkCommand::ProgressEvent, this, &vtkPVGeometryFilter::HandleGeometryFilterProgress);

  this->Controller = nullptr;
  this->SetController(vtkMultiProcessController::GetGlobalController());
  this->GenerateProcessIds = (this->Controller && this->Controller->GetNumberOfProcesses() > 1);

  this->OutlineSource = vtkOutlineSource::New();

  this->PassThroughCellIds = 1;
  this->PassThroughPointIds = 1;

  this->HideInternalAMRFaces = true;
  this->UseNonOverlappingAMRMetaDataForOutlines = true;
}

//----------------------------------------------------------------------------
vtkPVGeometryFilter::~vtkPVGeometryFilter()
{
  // Be careful how you delete these so that you don't foul up the garbage
  // collector.
  if (this->GeometryFilter)
  {
    vtkGeometryFilter* tmp = this->GeometryFilter;
    this->GeometryFilter = nullptr;
    tmp->Delete();
  }
  if (this->GenericGeometryFilter)
  {
    vtkGenericGeometryFilter* tmp = this->GenericGeometryFilter;
    this->GenericGeometryFilter = nullptr;
    tmp->Delete();
  }
  if (this->UnstructuredGridGeometryFilter)
  {
    vtkUnstructuredGridGeometryFilter* tmp = this->UnstructuredGridGeometryFilter;
    this->UnstructuredGridGeometryFilter = nullptr;
    tmp->Delete();
  }
  if (this->RecoverWireframeFilter)
  {
    vtkRecoverGeometryWireframe* tmp = this->RecoverWireframeFilter;
    this->RecoverWireframeFilter = nullptr;
    tmp->Delete();
  }
  if (this->FeatureEdgesFilter)
  {
    this->FeatureEdgesFilter->Delete();
  }
  this->OutlineSource->Delete();
  this->SetController(nullptr);
}

//----------------------------------------------------------------------------
int vtkPVGeometryFilter::RequestDataObject(
  vtkInformation*, vtkInformationVector** inputVector, vtkInformationVector* outputVector)
{
  vtkDataObject* input = vtkDataObject::GetData(inputVector[0], 0);
  vtkDataObject* output = vtkDataSet::GetData(outputVector, 0);

  if (input)
  {
    // If input is composite-data, then output is multi-block of polydata,
    // otherwise it's a poly data.
    if (vtkCompositeDataSet::SafeDownCast(input))
    {
      if (vtkMultiBlockDataSet::SafeDownCast(output) == nullptr)
      {
        if (vtkMultiBlockDataSet::SafeDownCast(input))
        {
          // Some developers have sub-classed vtkMultiBlockDataSet, in which
          // case, we try to preserve the type.
          output = input->NewInstance();
        }
        else
        {
          output = vtkMultiBlockDataSet::New();
        }
        outputVector->GetInformationObject(0)->Set(vtkDataObject::DATA_OBJECT(), output);
        output->FastDelete();
      }
      return 1;
    }

    if (vtkPolyData::SafeDownCast(output) == nullptr)
    {
      output = vtkPolyData::New();
      outputVector->GetInformationObject(0)->Set(vtkDataObject::DATA_OBJECT(), output);
      output->FastDelete();
    }
    return 1;
  }

  return 0;
}

//----------------------------------------------------------------------------
vtkExecutive* vtkPVGeometryFilter::CreateDefaultExecutive()
{
  return vtkCompositeDataPipeline::New();
}

//----------------------------------------------------------------------------
void vtkPVGeometryFilter::HandleGeometryFilterProgress(vtkObject* caller, unsigned long, void*)
{
  vtkAlgorithm* algorithm = vtkAlgorithm::SafeDownCast(caller);
  // This limits progress for only the GeometryFilter.
  if (algorithm)
  {
    double progress = algorithm->GetProgress();
    if (progress > 0.0 && progress < 1.0)
    {
      this->UpdateProgress(progress);
    }
    if (this->AbortExecute)
    {
      algorithm->SetAbortExecute(1);
    }
  }
}

//----------------------------------------------------------------------------
int vtkPVGeometryFilter::CheckAttributes(vtkDataObject* input)
{
  if (input->IsA("vtkDataSet"))
  {
    if (static_cast<vtkDataSet*>(input)->CheckAttributes())
    {
      return 1;
    }
  }
  else if (input->IsA("vtkCompositeDataSet"))
  {
    vtkCompositeDataSet* compInput = static_cast<vtkCompositeDataSet*>(input);
    vtkCompositeDataIterator* iter = compInput->NewIterator();
    iter->GoToFirstItem();
    while (!iter->IsDoneWithTraversal())
    {
      vtkDataObject* curDataSet = iter->GetCurrentDataObject();
      if (curDataSet && this->CheckAttributes(curDataSet))
      {
        return 1;
      }
      iter->GoToNextItem();
    }
    iter->Delete();
  }
  return 0;
}

//----------------------------------------------------------------------------
int vtkPVGeometryFilter::RequestUpdateExtent(
  vtkInformation* request, vtkInformationVector** inputVector, vtkInformationVector* outputVector)
{
  return this->Superclass::RequestUpdateExtent(request, inputVector, outputVector);
}

//----------------------------------------------------------------------------
void vtkPVGeometryFilter::ExecuteAMRBlockOutline(
  const double bounds[6], vtkPolyData* output, const bool extractface[6])
{
  // we generate outline faces, so that front face/back face culling works if
  // needed BUG #0011065. We only do this for AMR datasets for now, but we can
  // extend to all types of datasets, if needed.

  vtkNew<vtkPoints> points;
  points->Allocate(8);

  vtkNew<vtkCellArray> lines;
  lines->Allocate(lines->EstimateSize(12, 2));

  double x[3];
  x[0] = bounds[0];
  x[1] = bounds[2];
  x[2] = bounds[4];
  points->InsertPoint(0, x);
  x[0] = bounds[1];
  x[1] = bounds[2];
  x[2] = bounds[4];
  points->InsertPoint(1, x);
  x[0] = bounds[0];
  x[1] = bounds[3];
  x[2] = bounds[4];
  points->InsertPoint(2, x);
  x[0] = bounds[1];
  x[1] = bounds[3];
  x[2] = bounds[4];
  points->InsertPoint(3, x);
  x[0] = bounds[0];
  x[1] = bounds[2];
  x[2] = bounds[5];
  points->InsertPoint(4, x);
  x[0] = bounds[1];
  x[1] = bounds[2];
  x[2] = bounds[5];
  points->InsertPoint(5, x);
  x[0] = bounds[0];
  x[1] = bounds[3];
  x[2] = bounds[5];
  points->InsertPoint(6, x);
  x[0] = bounds[1];
  x[1] = bounds[3];
  x[2] = bounds[5];
  points->InsertPoint(7, x);

  // xmin face
  if (extractface[0])
  {
    vtkIdType pts[4] = { 0, 4, 6, 2 };
    lines->InsertNextCell(4, pts);
  }
  // xmax face
  if (extractface[1])
  {
    vtkIdType pts[4] = { 1, 3, 7, 5 };
    lines->InsertNextCell(4, pts);
  }

  // ymin face
  if (extractface[2])
  {
    vtkIdType pts[4] = { 0, 1, 5, 4 };
    lines->InsertNextCell(4, pts);
  }
  // ymax face
  if (extractface[3])
  {
    vtkIdType pts[4] = { 2, 6, 7, 3 };
    lines->InsertNextCell(4, pts);
  }

  // zmin face
  if (extractface[4])
  {
    vtkIdType pts[4] = { 0, 2, 3, 1 };
    lines->InsertNextCell(4, pts);
  }

  // zmax face
  if (extractface[5])
  {
    vtkIdType pts[4] = { 4, 5, 7, 6 };
    lines->InsertNextCell(4, pts);
  }

  output->SetPoints(points.GetPointer());
  output->SetPolys(lines.GetPointer());

  this->OutlineFlag = 1;
}

//----------------------------------------------------------------------------
void vtkPVGeometryFilter::ExecuteAMRBlock(
  vtkUniformGrid* input, vtkPolyData* output, const bool extractface[6])
{
  assert(input != nullptr && output != nullptr && this->UseOutline == 0);
  if (input->GetNumberOfCells() > 0)
  {
    int extent[6];
    input->GetExtent(extent);
    this->GeometryFilter->StructuredExecute(input, output, extent, const_cast<bool*>(extractface));
  }
  this->OutlineFlag = 0;
}

//----------------------------------------------------------------------------
void vtkPVGeometryFilter::ExecuteBlock(vtkDataObject* input, vtkPolyData* output, int doCommunicate,
  int updatePiece, int updateNumPieces, int updateGhosts, const int* wholeExtent)
{
  // Copy field data from the input block to the output block
  output->GetFieldData()->PassData(input->GetFieldData());

  if (auto imageData = vtkImageData::SafeDownCast(input))
  {
    this->ImageDataExecute(imageData, output, doCommunicate, updatePiece, wholeExtent);
  }
  else if (auto structuredGrid = vtkStructuredGrid::SafeDownCast(input))
  {
    this->StructuredGridExecute(
      structuredGrid, output, updatePiece, updateNumPieces, updateGhosts, wholeExtent);
  }
  else if (auto rectilinearGrid = vtkRectilinearGrid::SafeDownCast(input))
  {
    this->RectilinearGridExecute(
      rectilinearGrid, output, updatePiece, updateNumPieces, updateGhosts, wholeExtent);
  }
  else if (auto unstructuredGridBase = vtkUnstructuredGridBase::SafeDownCast(input))
  {
    this->UnstructuredGridExecute(unstructuredGridBase, output, doCommunicate);
  }
  else if (auto polyData = vtkPolyData::SafeDownCast(input))
  {
    this->PolyDataExecute(polyData, output, doCommunicate);
  }
  else if (auto hyperTreeGrid = vtkHyperTreeGrid::SafeDownCast(input))
  {
    this->HyperTreeGridExecute(hyperTreeGrid, output, doCommunicate);
  }
  else if (auto explicitStructuredGrid = vtkExplicitStructuredGrid::SafeDownCast(input))
  {
    this->ExplicitStructuredGridExecute(explicitStructuredGrid, output, doCommunicate, wholeExtent);
  }
  else if (auto dataset = vtkDataSet::SafeDownCast(input))
  {
    this->DataSetExecute(dataset, output, doCommunicate);
  }
  else if (auto genericDataSet = vtkGenericDataSet::SafeDownCast(input))
  {
    this->GenericDataSetExecute(genericDataSet, output, doCommunicate);
  }
  else if (auto cellGrid = vtkCellGrid::SafeDownCast(input))
  {
    this->CellGridExecute(cellGrid, output, doCommunicate);
  }
}

//----------------------------------------------------------------------------
int vtkPVGeometryFilter::RequestData(
  vtkInformation* request, vtkInformationVector** inputVector, vtkInformationVector* outputVector)
{
  this->GeometryFilter->SetRemoveGhostInterfaces(!this->GenerateFeatureEdges);
  vtkDataObject* input = vtkDataObject::GetData(inputVector[0], 0);
  if (vtkCompositeDataSet::SafeDownCast(input))
  {
    vtkTimerLog::MarkStartEvent("vtkPVGeometryFilter::RequestData");
    vtkGarbageCollector::DeferredCollectionPush();
    if (input->IsA("vtkUniformGridAMR"))
    {
      this->RequestAMRData(request, inputVector, outputVector);
    }
    else
    {
      this->RequestDataObjectTree(request, inputVector, outputVector);
    }
    vtkTimerLog::MarkStartEvent("vtkPVGeometryFilter::GarbageCollect");
    vtkGarbageCollector::DeferredCollectionPop();
    vtkTimerLog::MarkEndEvent("vtkPVGeometryFilter::GarbageCollect");
    vtkTimerLog::MarkEndEvent("vtkPVGeometryFilter::RequestData");
    return 1;
  }

  vtkPolyData* output = vtkPolyData::GetData(outputVector, 0);
  assert(output != nullptr);

  int procid = 0;
  int numProcs = 1;
  if (this->Controller)
  {
    procid = this->Controller->GetLocalProcessId();
    numProcs = this->Controller->GetNumberOfProcesses();
  }
  int* wholeExtent =
    vtkStreamingDemandDrivenPipeline::GetWholeExtent(inputVector[0]->GetInformationObject(0));
  this->ExecuteBlock(input, output, 1, procid, numProcs, 0, wholeExtent);
  this->CleanupOutputData(output, 1);
  return 1;
}

//----------------------------------------------------------------------------
void vtkPVGeometryFilter::CleanupOutputData(vtkPolyData* output, int doCommunicate)
{
  if (this->GenerateFeatureEdges)
  {
    this->FeatureEdgesFilter->SetInputData(output);
    this->FeatureEdgesFilter->Update();
    output->ShallowCopy(this->FeatureEdgesFilter->GetOutput());
  }
  this->ExecuteCellNormals(output, doCommunicate);
  output->RemoveGhostCells();
  if (this->GenerateProcessIds && output)
  {
    // add process ids array.
    unsigned int procId =
      this->Controller ? static_cast<unsigned int>(this->Controller->GetLocalProcessId()) : 0;

    vtkIdType numPoints = output->GetNumberOfPoints();
    if (numPoints > 0)
    {
      vtkNew<vtkUnsignedIntArray> array;
      array->SetNumberOfTuples(numPoints);
      array->FillTypedComponent(0, procId);
      array->SetName("vtkProcessId");
      output->GetPointData()->AddArray(array);
    }

    vtkIdType numCells = output->GetNumberOfCells();
    if (numCells > 0)
    {
      vtkNew<vtkUnsignedIntArray> cellArray;
      cellArray->SetNumberOfTuples(numCells);
      cellArray->FillTypedComponent(0, procId);
      cellArray->SetName("vtkProcessId");
      output->GetCellData()->AddArray(cellArray);
    }
  }
}

//----------------------------------------------------------------------------
namespace
{
static vtkPolyData* vtkPVGeometryFilterMergePieces(vtkPartitionedDataSet* mp)
{
  unsigned int num_pieces = mp->GetNumberOfPartitions();
  if (num_pieces == 0)
  {
    return nullptr;
  }

  std::vector<vtkPolyData*> inputs;
  std::vector<int> points_counts, cell_counts, verts_counts, polys_counts, lines_counts,
    strips_counts;

  polys_counts.resize(num_pieces);
  verts_counts.resize(num_pieces);
  lines_counts.resize(num_pieces);
  strips_counts.resize(num_pieces);
  points_counts.resize(num_pieces);
  cell_counts.resize(num_pieces);
  for (unsigned int cc = 0; cc < num_pieces; cc++)
  {
    vtkPolyData* piece = vtkPolyData::SafeDownCast(mp->GetPartition(cc));
    if (piece && piece->GetNumberOfPoints() > 0)
    {
      inputs.push_back(piece);
      points_counts[cc] = piece->GetNumberOfPoints();
      cell_counts[cc] = piece->GetNumberOfCells();
      verts_counts[cc] = piece->GetNumberOfVerts();
      polys_counts[cc] = piece->GetNumberOfPolys();
      lines_counts[cc] = piece->GetNumberOfLines();
      strips_counts[cc] = piece->GetNumberOfStrips();
    }
  }

  if (inputs.empty())
  {
    // not much to do, this is an empty multi-piece.
    return nullptr;
  }

  // Save field data attached to pieces as we want to preserve it (the vtkAppendPolydata filter
  // removes it). The field data is expected to be the same for all pieces of the multipiece so we
  // take it from the first one. See paraview/paraview#20748
  vtkFieldData* fd = inputs[0]->GetFieldData();

  vtkNew<vtkPolyData> output;
  vtkNew<vtkAppendPolyData> appender;
  appender->ExecuteAppend(output, inputs.data(), static_cast<int>(inputs.size()));
  inputs.clear();

  output->SetFieldData(fd);

  std::vector<int> points_offsets, verts_offsets, lines_offsets, polys_offsets, strips_offsets;
  polys_offsets.resize(num_pieces);
  verts_offsets.resize(num_pieces);
  lines_offsets.resize(num_pieces);
  strips_offsets.resize(num_pieces);
  points_offsets.resize(num_pieces);
  points_offsets[0] = 0;
  verts_offsets[0] = 0;
  lines_offsets[0] = output->GetNumberOfVerts();
  polys_offsets[0] = lines_offsets[0] + output->GetNumberOfLines();
  strips_offsets[0] = polys_offsets[0] + output->GetNumberOfPolys();
  for (unsigned int cc = 1; cc < num_pieces; cc++)
  {
    points_offsets[cc] = points_offsets[cc - 1] + points_counts[cc - 1];
    verts_offsets[cc] = verts_offsets[cc - 1] + verts_counts[cc - 1];
    lines_offsets[cc] = lines_offsets[cc - 1] + lines_counts[cc - 1];
    polys_offsets[cc] = polys_offsets[cc - 1] + polys_counts[cc - 1];
    strips_offsets[cc] = strips_offsets[cc - 1] + strips_counts[cc - 1];
  }

  for (unsigned int cc = 0; cc < num_pieces; cc++)
  {
    mp->SetPartition(cc, nullptr);
  }

  mp->SetPartition(0, output);

  vtkInformation* metadata = mp->GetMetaData(static_cast<unsigned int>(0));
  metadata->Set(
    vtkPVGeometryFilter::POINT_OFFSETS(), points_offsets.data(), static_cast<int>(num_pieces));
  metadata->Set(
    vtkPVGeometryFilter::VERTS_OFFSETS(), verts_offsets.data(), static_cast<int>(num_pieces));
  metadata->Set(
    vtkPVGeometryFilter::LINES_OFFSETS(), lines_offsets.data(), static_cast<int>(num_pieces));
  metadata->Set(
    vtkPVGeometryFilter::POLYS_OFFSETS(), polys_offsets.data(), static_cast<int>(num_pieces));
  metadata->Set(
    vtkPVGeometryFilter::STRIPS_OFFSETS(), strips_offsets.data(), static_cast<int>(num_pieces));
  return output;
}
};

//----------------------------------------------------------------------------
void vtkPVGeometryFilter::AddCompositeIndex(vtkPolyData* pd, unsigned int index)
{
  vtkNew<vtkUnsignedIntArray> cindex;
  cindex->SetNumberOfComponents(1);
  cindex->SetNumberOfTuples(pd->GetNumberOfCells());
  cindex->FillTypedComponent(0, index);
  cindex->SetName("vtkCompositeIndex");
  pd->GetCellData()->AddArray(cindex);

  vtkNew<vtkUnsignedIntArray> pindex;
  pindex->SetNumberOfComponents(1);
  pindex->SetNumberOfTuples(pd->GetNumberOfPoints());
  pindex->FillTypedComponent(0, index);
  pindex->SetName("vtkCompositeIndex");
  pd->GetPointData()->AddArray(pindex);
}

//----------------------------------------------------------------------------
void vtkPVGeometryFilter::AddBlockColors(vtkDataObject* pd, unsigned int index)
{
  vtkNew<vtkUnsignedIntArray> cindex;
  cindex->SetNumberOfComponents(1);
  cindex->SetNumberOfTuples(1);
  cindex->SetValue(0, index % this->BlockColorsDistinctValues);
  cindex->SetName("vtkBlockColors");
  pd->GetFieldData()->AddArray(cindex);
}

//----------------------------------------------------------------------------
void vtkPVGeometryFilter::AddHierarchicalIndex(
  vtkPolyData* pd, unsigned int level, unsigned int index)
{
  vtkNew<vtkUnsignedIntArray> dslevel;
  dslevel->SetNumberOfTuples(pd->GetNumberOfCells());
  dslevel->FillTypedComponent(0, level);
  dslevel->SetName("vtkAMRLevel");
  pd->GetCellData()->AddArray(dslevel);

  vtkNew<vtkUnsignedIntArray> dsindex;
  dsindex->SetNumberOfTuples(pd->GetNumberOfCells());
  dsindex->FillTypedComponent(0, index);
  dsindex->SetName("vtkAMRIndex");
  pd->GetCellData()->AddArray(dsindex);
}

//----------------------------------------------------------------------------
int vtkPVGeometryFilter::RequestAMRData(
  vtkInformation*, vtkInformationVector** inputVector, vtkInformationVector* outputVector)
{
  vtkTimerLog::MarkStartEvent("vtkPVGeometryFilter::RequestAMRData");

  // STEP 0: Acquire input & output object
  vtkMultiBlockDataSet* output = vtkMultiBlockDataSet::GetData(outputVector, 0);
  if (!output)
  {
    vtkErrorMacro("Output AMR multi-block dataset is NULL");
    return 0;
  }

  vtkUniformGridAMR* amr = vtkUniformGridAMR::GetData(inputVector[0], 0);
  if (!amr)
  {
    vtkErrorMacro("Input AMR dataset is NULL");
    return 0;
  }

  // STEP 1: Construct output object this will be multipiece that has all the
  // datasets under it.
  vtkNew<vtkMultiPieceDataSet> amrDatasets;
  output->SetNumberOfBlocks(1);
  output->SetBlock(0, amrDatasets.GetPointer());
  amrDatasets->SetNumberOfPieces(amr->GetTotalNumberOfBlocks());

  // STEP 2: Check Attributes
  vtkTimerLog::MarkStartEvent("vtkPVGeometryFilter::CheckAttributes");
  if (this->CheckAttributes(amr))
  {
    vtkErrorMacro("CheckAttributes() failed!");
    return 0;
  }
  vtkTimerLog::MarkEndEvent("vtkPVGeometryFilter::CheckAttributes");

  // STEP 3: Loop through data, determine if they are visible and call
  // execute block to get the polydata to render.
  // We use different implementations for overlapping and non-overlapping amrs.

  vtkOverlappingAMR* overlappingAMR = vtkOverlappingAMR::SafeDownCast(amr);

  double bounds[6];
  amr->GetBounds(bounds);
  if (this->Controller && this->Controller->GetNumberOfProcesses() > 1)
  {
    // Since bounds are not necessary synced up, especially for non-overlapping
    // AMR datasets, we sync them up across all processes.
    vtkPVGeometryFilter::BoundsReductionOperation operation;
    double received_bounds[6];
    this->Controller->AllReduce(bounds, received_bounds, 6, &operation);
    memcpy(bounds, received_bounds, sizeof(double) * 6);
  }

  unsigned int block_id = 0;
  for (unsigned int level = 0; level < amr->GetNumberOfLevels(); ++level)
  {
    unsigned int num_datasets = amr->GetNumberOfDataSets(level);
    for (unsigned int dataIdx = 0; dataIdx < num_datasets; ++dataIdx, block_id++)
    {
      vtkUniformGrid* ug = amr->GetDataSet(level, dataIdx);
      if ((ug == nullptr && this->UseOutline == 0) ||
        (ug == nullptr && this->UseOutline == 1 && overlappingAMR == nullptr))
      {
        // if this->UseOutline == 0,we need uniform grid to be present.

        // if this->UseOutline ==1, we need ug only for non-overlapping AMR. For
        // overlapping AMR, we can generate outline using the meta-data
        // available.
        continue;
      }

      if (overlappingAMR != nullptr && !this->UseNonOverlappingAMRMetaDataForOutlines &&
        ug == nullptr)
      {
        // for non-overlapping AMR, if we were told to not use meta-data, don't.
        continue;
      }

      double data_bounds[6];
      double error_margin = 0.01;

      // we have different mechanisms for determining if any of the faces of the
      // block are visible and what faces are visible based on the type of amr.
      if (overlappingAMR)
      {
        // for overlappingAMR, we use the meta-data to determine AMR bounds.
        overlappingAMR->GetAMRInfo()->GetBounds(level, dataIdx, data_bounds);
        double data_spacing[3];
        overlappingAMR->GetAMRInfo()->GetSpacing(level, data_spacing);
        error_margin = vtkMath::Norm(data_spacing);
      }
      else if (ug)
      {
        // for non-overlapping AMR, we use the bounds from the heavy-data
        // itself.
        ug->GetBounds(data_bounds);

        double data_spacing[3];
        ug->GetSpacing(data_spacing);
        error_margin = vtkMath::Norm(data_spacing);
      }
      else
      {
        continue; // skip block.
      }

      bool extractface[6] = { true, true, true, true, true, true };
      for (int cc = 0; this->HideInternalAMRFaces && cc < 6; cc++)
      {
        double delta = fabs(data_bounds[cc] - bounds[cc]);
        extractface[cc] = (delta < error_margin);
      }

      if (!(extractface[0] || extractface[1] || extractface[2] || extractface[3] ||
            extractface[4] || extractface[5]))
      {
        // we are not extracting a single face. nothing to do here.
        continue;
      }

      vtkNew<vtkPolyData> outputBlock;
      if (this->UseOutline)
      {
        this->ExecuteAMRBlockOutline(data_bounds, outputBlock.GetPointer(), extractface);
      }
      else
      {
        this->ExecuteAMRBlock(ug, outputBlock.GetPointer(), extractface);
      }
      if (!this->UseOutline)
      {
        // don't process attribute arrays when generating outlines.

        this->CleanupOutputData(outputBlock.GetPointer(), /*doCommunicate=*/0);
        this->AddCompositeIndex(outputBlock.GetPointer(), amr->GetCompositeIndex(level, dataIdx));
        this->AddHierarchicalIndex(outputBlock.GetPointer(), level, dataIdx);
        // we don't call this->AddBlockColors() for AMR dataset since it doesn't
        // make sense,  nor can be supported since all datasets merged into a
        // single polydata for rendering.
      }
      amrDatasets->SetPiece(block_id, outputBlock.GetPointer());
    }
  }

  // to avoid overburdening the rendering code with having to render a large
  // number of pieces, we merge the pieces.
  vtkPVGeometryFilterMergePieces(amrDatasets.GetPointer());

  // since we no longer care about the structure of the blocks in the composite
  // dataset (we are passing composite ids in the data itself to help identify
  // what block it came from), we can shrink allocated empty pointers for pieces
  // that vtkPVGeometryFilterMergePieces merged into one.
  amrDatasets->SetNumberOfPieces(1);

  vtkTimerLog::MarkEndEvent("vtkPVGeometryFilter::RequestAMRData");
  return 1;
}

//----------------------------------------------------------------------------
int vtkPVGeometryFilter::RequestDataObjectTree(
  vtkInformation*, vtkInformationVector** inputVector, vtkInformationVector* outputVector)
{
  vtkTimerLog::MarkStartEvent("vtkPVGeometryFilter::RequestDataObjectTree");

  vtkDataObjectTree* output = vtkDataObjectTree::GetData(outputVector, 0);
  if (!output)
  {
    return 0;
  }

  vtkDataObjectTree* input = vtkDataObjectTree::GetData(inputVector[0], 0);
  if (!input)
  {
    return 0;
  }
  output->CopyStructure(input);

  vtkTimerLog::MarkStartEvent("vtkPVGeometryFilter::CheckAttributes");
  if (this->CheckAttributes(input))
  {
    return 0;
  }
  vtkTimerLog::MarkEndEvent("vtkPVGeometryFilter::CheckAttributes");

  vtkTimerLog::MarkStartEvent("vtkPVGeometryFilter::ExecuteCompositeDataSet");
  vtkSmartPointer<vtkDataObjectTreeIterator> inIter;
  inIter.TakeReference(input->NewTreeIterator());
  inIter->VisitOnlyLeavesOn();
  inIter->SkipEmptyNodesOn();

  // get a block count for progress scaling.
  unsigned int totNumBlocks = 0;
  for (inIter->InitTraversal(); !inIter->IsDoneWithTraversal(); inIter->GoToNextItem())
  {
    ++totNumBlocks;
  }

  int* wholeExtent =
    vtkStreamingDemandDrivenPipeline::GetWholeExtent(inputVector[0]->GetInformationObject(0));
  int numInputs = 0;
  for (inIter->InitTraversal(); !inIter->IsDoneWithTraversal(); inIter->GoToNextItem())
  {
    vtkDataObject* block = inIter->GetCurrentDataObject();
    if (!block)
    {
      continue;
    }

    vtkNew<vtkPolyData> tmpOut;
    this->ExecuteBlock(block, tmpOut, 0, 0, 1, 0, wholeExtent);
    this->CleanupOutputData(tmpOut, 0);
    // skip empty nodes.
    if (tmpOut->GetNumberOfPoints() > 0)
    {
      output->SetDataSet(inIter, tmpOut);

      const unsigned int current_flat_index = inIter->GetCurrentFlatIndex();
      this->AddCompositeIndex(tmpOut, current_flat_index);
    }

    numInputs++;
    this->UpdateProgress(static_cast<float>(numInputs) / totNumBlocks);
  }
  vtkTimerLog::MarkEndEvent("vtkPVGeometryFilter::ExecuteCompositeDataSet");

  // Merge multi-pieces to avoid efficiency setbacks since multipieces can have
  // too many pieces.
  vtkSmartPointer<vtkDataObjectTreeIterator> outIter;
  outIter.TakeReference(output->NewTreeIterator());
  outIter->VisitOnlyLeavesOff();
  outIter->SkipEmptyNodesOn();

  std::vector<vtkPartitionedDataSet*> pieces_to_merge;
  for (outIter->InitTraversal(); !outIter->IsDoneWithTraversal(); outIter->GoToNextItem())
  {
    if (auto piece = vtkPartitionedDataSet::SafeDownCast(outIter->GetCurrentDataObject()))
    {
      pieces_to_merge.push_back(piece);
    }
  }

  // now merge these pieces (doing it in the above loop confuses the iterator).
  for (auto piece : pieces_to_merge)
  {
    vtkPVGeometryFilterMergePieces(piece);
  }

  if (this->Controller && this->Controller->GetNumberOfProcesses() > 1 &&
    !pieces_to_merge.empty() &&
    (vtkPartitionedDataSet::SafeDownCast(input) ||
      vtkPartitionedDataSetCollection::SafeDownCast(input)))
  {
    // since output of this filter is MB, and not a PDC or PD, we need to ensure
    // the number of pieces are same on all ranks (fixes #20654).
    std::vector<unsigned int> counts;
    counts.reserve(pieces_to_merge.size());
    for (auto& multipiece : pieces_to_merge)
    {
      counts.push_back(multipiece->GetNumberOfPartitions());
    }
    std::vector<unsigned int> gcounts(counts.size());
    this->Controller->AllReduce(counts.data(), gcounts.data(),
      static_cast<vtkIdType>(counts.size()), vtkCommunicator::MAX_OP);
    for (size_t cc = 0; cc < gcounts.size(); ++cc)
    {
      pieces_to_merge[cc]->SetNumberOfPartitions(gcounts[cc]);
    }
  }

  if (this->Controller && this->Controller->GetNumberOfProcesses() > 1)
  {
    // When running in parallel, processes may have nullptr-leaf nodes at
    // different locations. To make our life easier in subsequent filtering such as
    // vtkAllToNRedistributeCompositePolyData or vtkKdTreeManager we ensure that
    // all nullptr-leafs match up across processes i.e. if any leaf is non-nullptr on
    // any process, then all other processes add empty polydatas for that leaf.

    std::vector<unsigned char> non_null_leaves;
    non_null_leaves.reserve(totNumBlocks); // just an estimate.
    outIter->VisitOnlyLeavesOn();
    outIter->SkipEmptyNodesOn();
    for (outIter->InitTraversal(); !outIter->IsDoneWithTraversal(); outIter->GoToNextItem())
    {
      non_null_leaves.resize(outIter->GetCurrentFlatIndex() + 1, 0);
      non_null_leaves[outIter->GetCurrentFlatIndex()] = static_cast<unsigned char>(1);
    }

    int count = static_cast<int>(non_null_leaves.size());
    int reduced_size;
    this->Controller->AllReduce(&count, &reduced_size, 1, vtkCommunicator::MAX_OP);
    assert(reduced_size >= static_cast<int>(non_null_leaves.size()));
    non_null_leaves.resize(reduced_size, 0);
    // if reduced_size ==0, then all processes have no non-nullptr-leaves, so
    // nothing special to do here.
    if (reduced_size != 0)
    {
      std::vector<unsigned char> reduced_non_null_leaves;
      reduced_non_null_leaves.resize(reduced_size, 0);
      this->Controller->AllReduce(non_null_leaves.data(), reduced_non_null_leaves.data(),
        reduced_size, vtkCommunicator::MAX_OP);

      outIter->SkipEmptyNodesOff();
      outIter->VisitOnlyLeavesOn();
      for (outIter->InitTraversal(); !outIter->IsDoneWithTraversal(); outIter->GoToNextItem())
      {
        const unsigned int index = outIter->GetCurrentFlatIndex();
        if (outIter->GetCurrentDataObject() == nullptr &&
          index < static_cast<unsigned int>(reduced_non_null_leaves.size()) &&
          reduced_non_null_leaves[index] != 0)
        {
          vtkNew<vtkPolyData> trivalInput;
          this->AddCompositeIndex(trivalInput, index);
          output->SetDataSet(outIter, trivalInput);
        }
      }
    }
  }

  unsigned int block_id = 0;
  if (vtkPartitionedDataSetCollection::SafeDownCast(input))
  {
    // To avoid paraview/paraview#20908, we use a 2-level approach.
    vtkMultiBlockDataSet* outputMB = vtkMultiBlockDataSet::SafeDownCast(output);
    for (block_id = 0; block_id < outputMB->GetNumberOfBlocks(); ++block_id)
    {
      auto datasets = vtkCompositeDataSet::GetDataSets(outputMB->GetBlock(block_id));
      for (auto dataset : datasets)
      {
        this->AddBlockColors(dataset, block_id);
      }
    }
  }
  else
  {
    // At this point, all ranks have consistent tree structure with leaf nodes
    // non-nullptr at exactly same locations. This is a good point to assign block
    // colors.
    outIter->SkipEmptyNodesOff();
    outIter->VisitOnlyLeavesOn();
    for (outIter->InitTraversal(); !outIter->IsDoneWithTraversal();
         outIter->GoToNextItem(), ++block_id)
    {
      if (auto dobj = outIter->GetCurrentDataObject())
      {
        this->AddBlockColors(dobj, block_id);
      }
    }
  }

  if (block_id > 0)
  {
    // Add block colors to root-node's field data to keep it from being flagged as
    // partial.
    this->AddBlockColors(output, 0);
  }
  // vtkPVGeometryFilter does not support PDC as of now, therefore we need to encode the assembly
  // as a field data array in the output.
  if (auto inputPDC = vtkPartitionedDataSetCollection::SafeDownCast(input))
  {
    if (inputPDC->GetDataAssembly())
    {
      const auto dataAssemblyString = inputPDC->GetDataAssembly()->SerializeToXML(vtkIndent());
      vtkNew<vtkStringArray> dataAssemblyArray;
      dataAssemblyArray->SetName("vtkDataAssembly");
      dataAssemblyArray->InsertNextValue(dataAssemblyString.c_str());
      output->GetFieldData()->AddArray(dataAssemblyArray);
    }
  }

  vtkTimerLog::MarkEndEvent("vtkPVGeometryFilter::RequestDataObjectTree");
  return 1;
}

//----------------------------------------------------------------------------
// We need to change the mapper.  Now it always flat shades when cell normals
// are available.
void vtkPVGeometryFilter::ExecuteCellNormals(vtkPolyData* output, int doCommunicate)
{
  if (!this->GenerateCellNormals)
  {
    return;
  }

  // Do not generate cell normals if any of the processes
  // have lines, verts or strips.
  vtkCellArray* aPrim;
  int skip = 0;
  aPrim = output->GetVerts();
  if (aPrim && aPrim->GetNumberOfCells())
  {
    skip = 1;
  }
  aPrim = output->GetLines();
  if (aPrim && aPrim->GetNumberOfCells())
  {
    skip = 1;
  }
  aPrim = output->GetStrips();
  if (aPrim && aPrim->GetNumberOfCells())
  {
    skip = 1;
  }
  if (this->Controller && doCommunicate)
  {
    int reduced_skip = 0;
    if (!this->Controller->AllReduce(&skip, &reduced_skip, 1, vtkCommunicator::MAX_OP))
    {
      vtkErrorMacro("Failed to reduce correctly.");
      skip = 1;
    }
    else
    {
      skip = reduced_skip;
    }
  }
  if (skip)
  {
    return;
  }

  double polyNorm[3];
  vtkNew<vtkFloatArray> cellNormals;
  cellNormals->SetName("cellNormals");
  cellNormals->SetNumberOfComponents(3);
  cellNormals->Allocate(3 * output->GetNumberOfCells());

  aPrim = output->GetPolys();
  if (aPrim && aPrim->GetNumberOfCells())
  {
    vtkPoints* p = output->GetPoints();

    auto cellIter = vtk::TakeSmartPointer(aPrim->NewIterator());
    for (cellIter->GoToFirstCell(); !cellIter->IsDoneWithTraversal(); cellIter->GoToNextCell())
    {
      vtkIdList* cell = cellIter->GetCurrentCell();
      vtkPolygon::ComputeNormal(
        p, static_cast<int>(cell->GetNumberOfIds()), cell->GetPointer(0), polyNorm);
      cellNormals->InsertNextTuple(polyNorm);
    }
  }

  if (cellNormals->GetNumberOfTuples() != output->GetNumberOfCells())
  {
    vtkErrorMacro("Number of cell normals does not match output.");
    cellNormals->Delete();
    return;
  }

  output->GetCellData()->AddArray(cellNormals);
  output->GetCellData()->SetActiveNormals(cellNormals->GetName());
}

//----------------------------------------------------------------------------
void vtkPVGeometryFilter::DataSetExecute(vtkDataSet* input, vtkPolyData* output, int doCommunicate)
{
  double bds[6];
  int procid = 0;

  if (!doCommunicate && input->GetNumberOfPoints() == 0)
  {
    return;
  }

  if (this->Controller)
  {
    procid = this->Controller->GetLocalProcessId();
  }

  input->GetBounds(bds);

  vtkPVGeometryFilter::BoundsReductionOperation operation;
  if (procid && doCommunicate)
  {
    // Satellite node
    this->Controller->Reduce(bds, nullptr, 6, &operation, 0);
  }
  else
  {
    if (this->Controller && doCommunicate)
    {
      double tmp[6];
      this->Controller->Reduce(bds, tmp, 6, &operation, 0);
      memcpy(bds, tmp, 6 * sizeof(double));
    }

    if (bds[1] >= bds[0] && bds[3] >= bds[2] && bds[5] >= bds[4])
    {
      // only output in process 0.
      this->OutlineSource->SetBounds(bds);
      this->OutlineSource->Update();

      output->SetPoints(this->OutlineSource->GetOutput()->GetPoints());
      output->SetLines(this->OutlineSource->GetOutput()->GetLines());
    }
  }
}

//----------------------------------------------------------------------------
void vtkPVGeometryFilter::GenericDataSetExecute(
  vtkGenericDataSet* input, vtkPolyData* output, int doCommunicate)
{
  double bds[6];
  int procid = 0;

  if (!this->UseOutline)
  {
    this->OutlineFlag = 0;

    // Geometry filter
    this->GenericGeometryFilter->SetInputData(input);
    this->GenericGeometryFilter->Update();
    output->ShallowCopy(this->GenericGeometryFilter->GetOutput());

    return;
  }

  // Just outline
  this->OutlineFlag = 1;

  if (!doCommunicate && input->GetNumberOfPoints() == 0)
  {
    return;
  }

  if (this->Controller)
  {
    procid = this->Controller->GetLocalProcessId();
  }

  input->GetBounds(bds);

  vtkPVGeometryFilter::BoundsReductionOperation operation;
  if (procid && doCommunicate)
  {
    // Satellite node
    this->Controller->Reduce(bds, nullptr, 6, &operation, 0);
  }
  else
  {
    if (doCommunicate)
    {
      double tmp[6];
      this->Controller->Reduce(bds, tmp, 6, &operation, 0);
      memcpy(bds, tmp, 6 * sizeof(double));
    }

    // only output in process 0.
    this->OutlineSource->SetBounds(bds);
    this->OutlineSource->Update();

    output->SetPoints(this->OutlineSource->GetOutput()->GetPoints());
    output->SetLines(this->OutlineSource->GetOutput()->GetLines());
  }
}

//----------------------------------------------------------------------------
void vtkPVGeometryFilter::CellGridExecute(
  vtkCellGrid* input, vtkPolyData* output, int vtkNotUsed(doCommunicate))
{
  double bounds[6];
  input->GetBounds(bounds);
  vtkNew<vtkOutlineSource> outline;
  outline->SetBounds(bounds);
  outline->Update();

  output->SetPoints(outline->GetOutput()->GetPoints());
  output->SetLines(outline->GetOutput()->GetLines());
  output->SetPolys(outline->GetOutput()->GetPolys());
}

//----------------------------------------------------------------------------
void vtkPVGeometryFilter::ImageDataExecute(
  vtkImageData* input, vtkPolyData* output, int doCommunicate, int updatePiece, const int* ext)
{
  double* spacing;
  double* origin;
  //   int* ext;
  double bounds[6];

  // If doCommunicate is false, use extent because the block is
  // entirely contained in this process.
  if (!doCommunicate)
  {
    ext = input->GetExtent();
  }

  // If 2d then default to superclass behavior.
  //  if (ext[0] == ext[1] || ext[2] == ext[3] || ext[4] == ext[5] ||
  //      !this->UseOutline)
  if (!this->UseOutline)
  {
    if (input->GetNumberOfCells() > 0)
    {
      this->GeometryFilter->StructuredExecute(
        input, output, const_cast<int*>(ext), nullptr, nullptr);
    }
    this->OutlineFlag = 0;
    return;
  }
  this->OutlineFlag = 1;

  //
  // Otherwise, let OutlineSource do all the work
  //

  if (ext[1] >= ext[0] && ext[3] >= ext[2] && ext[5] >= ext[4] &&
    (updatePiece == 0 || !doCommunicate))
  {
    spacing = input->GetSpacing();
    origin = input->GetOrigin();

    bounds[0] = spacing[0] * ((float)ext[0]) + origin[0];
    bounds[1] = spacing[0] * ((float)ext[1]) + origin[0];
    bounds[2] = spacing[1] * ((float)ext[2]) + origin[1];
    bounds[3] = spacing[1] * ((float)ext[3]) + origin[1];
    bounds[4] = spacing[2] * ((float)ext[4]) + origin[2];
    bounds[5] = spacing[2] * ((float)ext[5]) + origin[2];

    vtkNew<vtkOutlineSource> outline;
    outline->SetBounds(bounds);
    outline->Update();

    output->SetPoints(outline->GetOutput()->GetPoints());
    output->SetLines(outline->GetOutput()->GetLines());
    output->SetPolys(outline->GetOutput()->GetPolys());
  }
  else
  {
    vtkNew<vtkPoints> pts;
    output->SetPoints(pts);
  }
}

//----------------------------------------------------------------------------
void vtkPVGeometryFilter::StructuredGridExecute(vtkStructuredGrid* input, vtkPolyData* output,
  int updatePiece, int updateNumPieces, int updateGhosts, const int* wholeExtentArg)
{
  int wholeExtent[6];
  ::GetValidWholeExtent(input, wholeExtentArg, wholeExtent);

  if (!this->UseOutline)
  {
    if (input->GetNumberOfCells() > 0)
    {
      this->GeometryFilter->StructuredExecute(input, output, wholeExtent, nullptr, nullptr);
    }
    this->OutlineFlag = 0;
    return;
  }
  this->OutlineFlag = 1;

  vtkNew<vtkPVTrivialProducer> producer;
  producer->SetOutput(input);
  producer->SetWholeExtent(
    wholeExtent[0], wholeExtent[1], wholeExtent[2], wholeExtent[3], wholeExtent[4], wholeExtent[5]);

  vtkNew<vtkStructuredGridOutlineFilter> outline;
  outline->SetInputConnection(producer->GetOutputPort());
  outline->UpdatePiece(updatePiece, updateNumPieces, updateGhosts);
  output->CopyStructure(outline->GetOutput());
}

//----------------------------------------------------------------------------
void vtkPVGeometryFilter::RectilinearGridExecute(vtkRectilinearGrid* input, vtkPolyData* output,
  int vtkNotUsed(updatePiece), int vtkNotUsed(updateNumPieces), int vtkNotUsed(updateGhosts),
  const int* wholeExtent)
{
  if (!this->UseOutline)
  {
    if (input->GetNumberOfCells() > 0)
    {
      this->GeometryFilter->StructuredExecute(
        input, output, const_cast<int*>(wholeExtent), nullptr, nullptr);
    }
    this->OutlineFlag = 0;
    return;
  }
  this->OutlineFlag = 1;

  vtkNew<vtkPVTrivialProducer> producer;
  producer->SetOutput(input);
  producer->SetWholeExtent(
    wholeExtent[0], wholeExtent[1], wholeExtent[2], wholeExtent[3], wholeExtent[4], wholeExtent[5]);

  vtkNew<vtkRectilinearGridOutlineFilter> outline;
  outline->SetInputConnection(producer->GetOutputPort());
  outline->Update();
  output->CopyStructure(outline->GetOutput());
}

//----------------------------------------------------------------------------
void vtkPVGeometryFilter::UnstructuredGridExecute(
  vtkUnstructuredGridBase* input, vtkPolyData* output, int doCommunicate)
{
  if (!this->UseOutline)
  {
    this->OutlineFlag = 0;

    bool handleSubdivision = (this->Triangulate != 0) && (input->GetNumberOfCells() > 0);
    if (!handleSubdivision && (this->NonlinearSubdivisionLevel > 0))
    {
      // Check to see if the data actually has nonlinear cells.  Handling
      // nonlinear cells adds unnecessary work if we only have linear cells.
      if (input->GetNumberOfCells() > 0)
      {
        auto helper = vtkGeometryFilterHelper::CharacterizeUnstructuredGrid(input);
        if (!helper->IsLinear)
        {
          handleSubdivision = true;
        }
        delete helper;
      }
    }

    vtkSmartPointer<vtkIdTypeArray> facePtIds2OriginalPtIds;

    auto inputClone = vtkSmartPointer<vtkUnstructuredGridBase>::Take(input->NewInstance());
    inputClone->ShallowCopy(input);
    input = inputClone;

    if (handleSubdivision)
    {
      // Use the vtkUnstructuredGridGeometryFilter to extract 2D surface cells
      // from the geometry.  This is important to extract an appropriate
      // wireframe in vtkRecoverGeometryWireframe.  Also, at the time of this
      // writing vtkGeometryFilter only properly subdivides 2D cells past level 1.
      this->UnstructuredGridGeometryFilter->SetInputData(input);

      // Let the vtkUnstructuredGridGeometryFilter record from which point and
      // cell each face comes from in the standard vtkOriginalCellIds array.
      this->UnstructuredGridGeometryFilter->SetPassThroughCellIds(this->PassThroughCellIds);
      this->UnstructuredGridGeometryFilter->SetPassThroughPointIds(this->PassThroughPointIds);

      this->UnstructuredGridGeometryFilter->SetMatchBoundariesIgnoringCellOrder(
        this->MatchBoundariesIgnoringCellOrder);

      // Turn off ghost cell clipping. This ensures that ghost cells are retained
      // and handed to the GeometryFilter to ensure only valid faces are
      // generated. If this weren't here, then the GeometryFilter would
      // generate boundary faces between normal cells and where the ghost cells
      // used to be, which is not correct.
      this->UnstructuredGridGeometryFilter->DuplicateGhostCellClippingOff();

      // Disable point merging as it may prevent the correct visualization
      // of non-continuous attributes.
      this->UnstructuredGridGeometryFilter->MergingOff();

      // TODO: Make the consecutive internal filter execution have monotonically
      // increasing progress rather than restarting for every internal filter.
      this->UnstructuredGridGeometryFilter->Update();

      this->UnstructuredGridGeometryFilter->SetInputData(nullptr);

      // Feed the extracted surface as the input to the rest of the processing.
      input->ShallowCopy(this->UnstructuredGridGeometryFilter->GetOutput());

      // Keep a handle to the vtkOriginalPointIds array.  We might need it.
      facePtIds2OriginalPtIds =
        vtkIdTypeArray::SafeDownCast(input->GetPointData()->GetArray("vtkOriginalPointIds"));

      // Flag the data set surface filter to record original cell ids, but do it
      // in a specially named array that vtkRecoverGeometryWireframe will later
      // use.  Note that because the data set comes from
      // UnstructuredGridGeometryFilter, the ids will represent the faces rather
      // than the original cells, which is important.
      this->GeometryFilter->PassThroughCellIdsOn();
      this->GeometryFilter->SetOriginalCellIdsName(details::ORIGINAL_FACE_IDS);

      if (this->PassThroughPointIds)
      {
        // vtkGeometryFilter is going to strip the vtkOriginalPointIds
        // created by the vtkPVUnstructuredGridGeometryFilter because it
        // cannot interpolate the ids.  Make the vtkGeometryFilter make
        // its own original ids array.  We will resolve them later.
        this->GeometryFilter->PassThroughPointIdsOn();
      }
    }

    if (input->GetNumberOfCells() > 0)
    {
      this->GeometryFilter->UnstructuredGridExecute(input, output);
    }

    if (this->Triangulate && (output->GetNumberOfPolys() > 0))
    {
      // Triangulate the polygonal mesh if requested to avoid rendering
      // issues of non-convex polygons.
      vtkNew<vtkTriangleFilter> triangleFilter;
      triangleFilter->SetInputData(output);
      triangleFilter->Update();
      output->ShallowCopy(triangleFilter->GetOutput());
    }

    if (handleSubdivision && !this->GenerateFeatureEdges)
    {
      // Restore state of GeometryFilter.
      this->GeometryFilter->SetPassThroughCellIds(this->PassThroughCellIds);
      this->GeometryFilter->SetOriginalCellIdsName(nullptr);
      this->GeometryFilter->SetPassThroughPointIds(this->PassThroughPointIds);

      this->GeometryFilter->SetMatchBoundariesIgnoringCellOrder(
        this->MatchBoundariesIgnoringCellOrder);

      // Now use vtkRecoverGeometryWireframe to create an edge flag attribute
      // that will cause the wireframe to be rendered correctly.
      vtkNew<vtkPolyData> nextStageInput;
      nextStageInput->ShallowCopy(output); // Yes output is correct.
      this->RecoverWireframeFilter->SetInputData(nextStageInput.Get());
      this->RecoverWireframeFilter->SetCellIdsAttribute(details::ORIGINAL_FACE_IDS);
      // TODO: Make the consecutive internal filter execution have monotonically
      // increasing progress rather than restarting for every internal filter.
      this->RecoverWireframeFilter->Update();
      this->RecoverWireframeFilter->SetInputData(nullptr);

      // Get what should be the final output.
      output->ShallowCopy(this->RecoverWireframeFilter->GetOutput());

      if (this->PassThroughPointIds)
      {
        // The output currently has a vtkOriginalPointIds array that maps points
        // to the data containing only the faces.  Correct this to point to the
        // original data set.
        vtkIdTypeArray* polyPtIds2FacePtIds =
          vtkIdTypeArray::SafeDownCast(output->GetPointData()->GetArray("vtkOriginalPointIds"));
        if (!facePtIds2OriginalPtIds || !polyPtIds2FacePtIds)
        {
          vtkErrorMacro(<< "Missing original point id arrays.");
          return;
        }
        vtkIdType numPts = polyPtIds2FacePtIds->GetNumberOfTuples();
        vtkNew<vtkIdTypeArray> polyPtIds2OriginalPtIds;
        polyPtIds2OriginalPtIds->SetName("vtkOriginalPointIds");
        polyPtIds2OriginalPtIds->SetNumberOfComponents(1);
        polyPtIds2OriginalPtIds->SetNumberOfTuples(numPts);
        for (vtkIdType polyPtId = 0; polyPtId < numPts; polyPtId++)
        {
          vtkIdType facePtId = polyPtIds2FacePtIds->GetValue(polyPtId);
          vtkIdType originalPtId = -1;
          if (facePtId >= 0)
          {
            originalPtId = facePtIds2OriginalPtIds->GetValue(facePtId);
          }
          polyPtIds2OriginalPtIds->SetValue(polyPtId, originalPtId);
        }
        output->GetPointData()->AddArray(polyPtIds2OriginalPtIds.Get());
      }
    }

    output->GetCellData()->RemoveArray(details::ORIGINAL_FACE_IDS);
    return;
  }

  this->OutlineFlag = 1;

  this->DataSetExecute(input, output, doCommunicate);
}

//----------------------------------------------------------------------------
void vtkPVGeometryFilter::PolyDataExecute(
  vtkPolyData* input, vtkPolyData* output, int doCommunicate)
{
  if (!this->UseOutline)
  {
    this->OutlineFlag = 0;
    output->ShallowCopy(input);
    if (this->PassThroughCellIds)
    {
      vtkNew<vtkIdTypeArray> originalCellIds;
      originalCellIds->SetName("vtkOriginalCellIds");
      originalCellIds->SetNumberOfComponents(1);
      vtkNew<vtkIdTypeArray> originalFaceIds;
      originalFaceIds->SetName(details::ORIGINAL_FACE_IDS);
      originalFaceIds->SetNumberOfComponents(1);
      vtkCellData* outputCD = output->GetCellData();
      outputCD->AddArray(originalCellIds.Get());
      if (this->Triangulate)
      {
        outputCD->AddArray(originalFaceIds.Get());
      }
      vtkIdType numTup = output->GetNumberOfCells();
      originalCellIds->SetNumberOfValues(numTup);
      originalFaceIds->SetNumberOfValues(numTup);
      for (vtkIdType cId = 0; cId < numTup; cId++)
      {
        originalCellIds->SetValue(cId, cId);
        originalFaceIds->SetValue(cId, cId);
      }
    }
    if (this->PassThroughPointIds)
    {
      vtkNew<vtkIdTypeArray> originalPointIds;
      originalPointIds->SetName("vtkOriginalPointIds");
      originalPointIds->SetNumberOfComponents(1);
      vtkPointData* outputPD = output->GetPointData();
      outputPD->AddArray(originalPointIds.Get());
      vtkIdType numTup = output->GetNumberOfPoints();
      originalPointIds->SetNumberOfValues(numTup);
      for (vtkIdType pId = 0; pId < numTup; pId++)
      {
        originalPointIds->SetValue(pId, pId);
      }
    }

    if (this->Triangulate)
    {
      // Triangulate the polygonal mesh.
      vtkNew<vtkTriangleFilter> triangleFilter;
      triangleFilter->SetInputData(output);
      triangleFilter->Update();

      // Now use vtkRecoverGeometryWireframe to create an edge flag attribute
      // that will cause the wireframe to be rendered correctly.
      this->RecoverWireframeFilter->SetInputData(triangleFilter->GetOutput());
      // TODO: Make the consecutive internal filter execution have monotonically
      // increasing progress rather than restarting for every internal filter.
      this->RecoverWireframeFilter->Update();
      this->RecoverWireframeFilter->SetInputData(nullptr);

      // Get what should be the final output.
      output->ShallowCopy(this->RecoverWireframeFilter->GetOutput());

      output->GetCellData()->RemoveArray(details::ORIGINAL_FACE_IDS);
    }
    return;
  }

  this->OutlineFlag = 1;
  this->DataSetExecute(input, output, doCommunicate);
}

//----------------------------------------------------------------------------
void vtkPVGeometryFilter::HyperTreeGridExecute(
  vtkHyperTreeGrid* input, vtkPolyData* output, int doCommunicate)
{
  if (!this->UseOutline)
  {
    this->OutlineFlag = 0;

    vtkNew<vtkHyperTreeGridGeometry> internalFilter;
    vtkNew<vtkHyperTreeGrid> htgCopy;
    htgCopy->ShallowCopy(input);
    internalFilter->SetInputData(htgCopy);
    internalFilter->SetPassThroughCellIds(this->PassThroughCellIds);
    internalFilter->SetOriginalCellIdArrayName("vtkOriginalCellIds");
    internalFilter->Update();
    output->ShallowCopy(internalFilter->GetOutput());
    return;
  }

  this->OutlineFlag = 1;
  double bds[6];
  int procid = 0;
  if (!doCommunicate && input->GetNumberOfCells() == 0)
  {
    return;
  }

  if (this->Controller)
  {
    procid = this->Controller->GetLocalProcessId();
  }

  input->GetBounds(bds);

  vtkPVGeometryFilter::BoundsReductionOperation operation;
  if (procid && doCommunicate)
  {
    // Satellite node
    this->Controller->Reduce(bds, nullptr, 6, &operation, 0);
  }
  else
  {
    if (this->Controller && doCommunicate)
    {
      double tmp[6];
      this->Controller->Reduce(bds, tmp, 6, &operation, 0);
      memcpy(bds, tmp, 6 * sizeof(double));
    }

    if (bds[1] >= bds[0] && bds[3] >= bds[2] && bds[5] >= bds[4])
    {
      // only output in process 0.
      this->OutlineSource->SetBounds(bds);
      this->OutlineSource->Update();

      output->SetPoints(this->OutlineSource->GetOutput()->GetPoints());
      output->SetLines(this->OutlineSource->GetOutput()->GetLines());
    }
  }
}

//----------------------------------------------------------------------------
void vtkPVGeometryFilter::ExplicitStructuredGridExecute(
  vtkExplicitStructuredGrid* input, vtkPolyData* out, int doCommunicate, const int* wholeExtent)
{
  vtkNew<vtkPVTrivialProducer> producer;
  producer->SetOutput(input);
  producer->SetWholeExtent(
    wholeExtent[0], wholeExtent[1], wholeExtent[2], wholeExtent[3], wholeExtent[4], wholeExtent[5]);
  producer->Update();

  if (!this->UseOutline)
  {
    this->OutlineFlag = 0;

    vtkNew<vtkExplicitStructuredGridSurfaceFilter> internalFilter;
    internalFilter->SetPassThroughPointIds(this->PassThroughPointIds);
    internalFilter->SetPassThroughCellIds(this->PassThroughCellIds);
    internalFilter->SetInputConnection(producer->GetOutputPort());
    internalFilter->Update();
    out->ShallowCopy(internalFilter->GetOutput());
    return;
  }
  vtkExplicitStructuredGrid* in =
    vtkExplicitStructuredGrid::SafeDownCast(producer->GetOutputDataObject(0));

  this->OutlineFlag = 1;
  this->DataSetExecute(in, out, doCommunicate);
}

//----------------------------------------------------------------------------
int vtkPVGeometryFilter::FillInputPortInformation(int port, vtkInformation* info)
{
  if (!this->Superclass::FillInputPortInformation(port, info))
  {
    return 0;
  }

  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
  info->Append(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkGenericDataSet");
  info->Append(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkCompositeDataSet");
  info->Append(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkHyperTreeGrid");
  info->Append(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkCellGrid");
  return 1;
}

//-----------------------------------------------------------------------------
void vtkPVGeometryFilter::ReportReferences(vtkGarbageCollector* collector)
{
  this->Superclass::ReportReferences(collector);
  vtkGarbageCollectorReport(collector, this->GeometryFilter, "GeometryFilter");
  vtkGarbageCollectorReport(collector, this->GenericGeometryFilter, "GenericGeometryFilter");
  vtkGarbageCollectorReport(
    collector, this->UnstructuredGridGeometryFilter, "UnstructuredGridGeometryFilter");
  vtkGarbageCollectorReport(collector, this->RecoverWireframeFilter, "RecoverWireframeFilter");
}

//----------------------------------------------------------------------------
void vtkPVGeometryFilter::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);

  os << indent << "OutlineFlag: " << (this->OutlineFlag ? "on" : "off") << endl;
  os << indent << "UseOutline: " << (this->UseOutline ? "on" : "off") << endl;
  os << indent << "GenerateFeatureEdges: " << (this->GenerateFeatureEdges ? "on" : "off") << endl;
  os << indent << "BlockColorsDistinctValues: " << this->BlockColorsDistinctValues << endl;
  os << indent << "GenerateCellNormals: " << (this->GenerateCellNormals ? "on" : "off") << endl;
  os << indent << "Triangulate: " << (this->Triangulate ? "on" : "off") << endl;
  os << indent << "NonlinearSubdivisionLevel: " << this->NonlinearSubdivisionLevel << endl;
  os << indent << "MatchBoundariesIgnoringCellOrder: " << this->MatchBoundariesIgnoringCellOrder
     << endl;
  os << indent << "Controller: " << this->Controller << endl;
  os << indent << "PassThroughCellIds: " << (this->PassThroughCellIds ? "on" : "off") << endl;
  os << indent << "PassThroughPointIds: " << (this->PassThroughPointIds ? "on" : "off") << endl;
  os << indent << "GenerateProcessIds: " << (this->GenerateProcessIds ? "on" : "off") << endl;
  os << indent << "HideInternalAMRFaces: " << (this->HideInternalAMRFaces ? "on" : "off") << endl;
  os << indent << "UseNonOverlappingAMRMetaDataForOutlines: "
     << (this->UseNonOverlappingAMRMetaDataForOutlines ? "on" : "off") << endl;
}

//----------------------------------------------------------------------------
void vtkPVGeometryFilter::SetPassThroughCellIds(int newvalue)
{
  this->PassThroughCellIds = newvalue;
  if (this->GeometryFilter)
  {
    this->GeometryFilter->SetPassThroughCellIds(this->PassThroughCellIds);
  }
  if (this->GenericGeometryFilter)
  {
    this->GenericGeometryFilter->SetPassThroughCellIds(this->PassThroughCellIds);
  }
}

//----------------------------------------------------------------------------
void vtkPVGeometryFilter::SetPassThroughPointIds(int newvalue)
{
  this->PassThroughPointIds = newvalue;
  if (this->GeometryFilter)
  {
    this->GeometryFilter->SetPassThroughPointIds(this->PassThroughPointIds);
  }
}

//-----------------------------------------------------------------------------
void vtkPVGeometryFilter::SetNonlinearSubdivisionLevel(int newvalue)
{
  if (this->NonlinearSubdivisionLevel != newvalue)
  {
    this->NonlinearSubdivisionLevel = newvalue;

    if (this->GeometryFilter)
    {
      this->GeometryFilter->SetNonlinearSubdivisionLevel(this->NonlinearSubdivisionLevel);
    }

    this->Modified();
  }
}

//-----------------------------------------------------------------------------
void vtkPVGeometryFilter::SetMatchBoundariesIgnoringCellOrder(int newvalue)
{
  if (this->MatchBoundariesIgnoringCellOrder != newvalue)
  {
    this->MatchBoundariesIgnoringCellOrder = newvalue;

    if (this->GeometryFilter)
    {
      this->GeometryFilter->SetMatchBoundariesIgnoringCellOrder(
        this->MatchBoundariesIgnoringCellOrder);
    }

    this->Modified();
  }
}
