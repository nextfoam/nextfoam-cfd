//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#ifndef vtk_m_filter_field_transform_PointElevation_h
#define vtk_m_filter_field_transform_PointElevation_h

#include <vtkm/filter/NewFilterField.h>
#include <vtkm/filter/field_transform/vtkm_filter_field_transform_export.h>

namespace vtkm
{
namespace filter
{
namespace field_transform
{
/// \brief  generate a scalar field along a specified direction
///
/// Generate scalar field from a dataset.
/// The scalar field values lie within a user specified range, and
/// are generated by computing a projection of each dataset point onto
/// a line. The line can be oriented arbitrarily. A typical example is
/// to generate scalars based on elevation or height above a plane.
class VTKM_FILTER_FIELD_TRANSFORM_EXPORT PointElevation : public vtkm::filter::NewFilterField
{
public:
  VTKM_CONT
  PointElevation();

  VTKM_CONT void SetLowPoint(const vtkm::Vec3f_64& point) { this->LowPoint = point; }
  VTKM_CONT void SetLowPoint(vtkm::Float64 x, vtkm::Float64 y, vtkm::Float64 z)
  {
    this->SetLowPoint({ x, y, z });
  }

  VTKM_CONT void SetHighPoint(const vtkm::Vec3f_64& point) { this->HighPoint = point; }
  VTKM_CONT void SetHighPoint(vtkm::Float64 x, vtkm::Float64 y, vtkm::Float64 z)
  {
    this->SetHighPoint({ x, y, z });
  }

  VTKM_CONT
  void SetRange(vtkm::Float64 low, vtkm::Float64 high)
  {
    this->RangeLow = low;
    this->RangeHigh = high;
  }

private:
  VTKM_CONT vtkm::cont::DataSet DoExecute(const vtkm::cont::DataSet& input) override;

  vtkm::Vec3f_64 LowPoint = { 0.0, 0.0, 0.0 };
  vtkm::Vec3f_64 HighPoint = { 0.0, 0.0, 1.0 };
  vtkm::Float64 RangeLow = 0.0, RangeHigh = 1.0;
};
} // namespace field_transform
class VTKM_DEPRECATED(1.8, "Use vtkm::filter::field_transform::PointElevation.") PointElevation
  : public vtkm::filter::field_transform::PointElevation
{
  using field_transform::PointElevation::PointElevation;
};
} // namespace filter
} // namespace vtkm

#endif // vtk_m_filter_field_transform_PointElevation_h