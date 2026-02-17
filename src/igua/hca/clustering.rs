use pyo3::prelude::*;
use pyo3::exceptions::PyValueError;
use pyo3::exceptions::PyTypeError;
use numpy::PyArray;
use numpy::Element;
use numpy::PyArrayMethods;
use numpy::PyUntypedArrayMethods;
use kodama::Method;
use num_traits::Float;
use num_traits::ToPrimitive;
use half::f16;


fn extract_dendrogram<'py, T>(
    py: Python<'py>, 
    dendrogram: kodama::Dendrogram<T>
) -> PyResult<Bound<'py, PyArray<f64, numpy::Ix2>>>
where
    T: Float,
{
    let steps = dendrogram.steps();
    unsafe {
        let z = PyArray::new(py, [steps.len(), 4], false);
        let z_view = z.try_readwrite()?;
        for (i, step) in dendrogram.steps().iter().enumerate() {
            *z_view.uget_mut([i, 0]) = step.cluster1.to_f64().unwrap();
            *z_view.uget_mut([i, 1]) = step.cluster2.to_f64().unwrap();
            *z_view.uget_mut([i, 2]) = step.dissimilarity.to_f64().unwrap();
            *z_view.uget_mut([i, 3]) = step.size.to_f64().unwrap();
        }
        Ok(z)
    }
}

fn linkage_impl<'py, T>(
    py: Python<'py>,
    distances: &Bound<'py, PyArray<T, numpy::Ix1>>,
    method: kodama::Method,
) -> PyResult<Bound<'py, PyArray<f64, numpy::Ix2>>>
where
    T: Element + Float,
{
    let n = (distances.len() as f32 * 2.0).sqrt().ceil() as usize;
    let mut readwrite = distances.try_readwrite()?;
    let slice = readwrite.as_slice_mut()?;
    let dendrogram = kodama::linkage(slice, n, method);
    extract_dendrogram(py, dendrogram)
}

#[pyfunction]
#[pyo3(signature = (distances, method="single"))]
pub fn linkage<'py>(
    py: Python<'py>,
    distances: &Bound<'py, PyAny>,
    method: &str,
) -> PyResult<Bound<'py, PyArray<f64, numpy::Ix2>>> {
    let variant = match method {
        "single" => Method::Single,
        "complete" => Method::Complete,
        "average" => Method::Average,
        "weighted" => Method::Weighted,
        "ward" => Method::Ward,
        "centroid" => Method::Centroid,
        "median" => Method::Median,
        other => return Err(PyValueError::new_err(format!("Invalid method: {}", other))),
    };

    let d = distances.as_borrowed();
    if let Ok(d) = <Bound<'py, PyArray::<f64, numpy::Ix1>>>::extract(d) {
        return linkage_impl(py, &d, variant);
    }
    if let Ok(d) = <Bound<'py, PyArray::<f32, numpy::Ix1>>>::extract(d) {
        return linkage_impl(py, &d, variant);
    }
    if let Ok(d) = <Bound<'py, PyArray::<f16, numpy::Ix1>>>::extract(d) {
        return linkage_impl(py, &d, variant);
    }

    Err(PyTypeError::new_err("Unsupported dtype in `linkage`"))
}
