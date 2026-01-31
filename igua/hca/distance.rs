use std::cmp::Ord;
use std::cmp::Ordering;

use pyo3::prelude::*;
use pyo3::exceptions::PyTypeError;
use rayon::prelude::*;
use half::f16;
use numpy::Element;
use numpy::PyArray;
use numpy::PyArrayMethods;
use num_traits::Float;
use num_traits::FromPrimitive;

fn manhattan_pair_impl<T>(
    indptr: &[i32],
    indices: &[i32],
    data: &[i32],
    weights: Option<&[i32]>,
    px: usize,
    py: usize,
) -> T 
where
    T: Float + Element + FromPrimitive,
{
    let i_next = indptr[px + 1] as usize;
    let j_next = indptr[py + 1] as usize;

    let mut i = indptr[px] as usize;
    let mut j = indptr[py] as usize;
    let mut d = 0;

    while i < i_next && j < j_next {
        match indices[i].cmp(&indices[j]) {
            Ordering::Equal => {
                let w = weights.map(|w| w[indices[i] as usize]).unwrap_or(1);
                d += data[i].abs_diff(data[j]) as u64 * w as u64;
                i += 1;
                j += 1;
            }
            Ordering::Less => {
                let w = weights.map(|w| w[indices[i] as usize]).unwrap_or(1);
                d += data[i].abs() as u64 * w as u64;
                i += 1;
            }
            Ordering::Greater => {
                let w = weights.map(|w| w[indices[j] as usize]).unwrap_or(1);
                d += data[j].abs() as u64 * w as u64;
                j += 1;
            }
        }
    }

    if i == i_next {
        while j < j_next {
                let w = weights.map(|w| w[indices[j] as usize]).unwrap_or(1);
            d += data[j].abs() as u64 * w as u64;
            j += 1;
        }
    } else {
        while i < i_next {
            let w = weights.map(|w| w[indices[i] as usize]).unwrap_or(1);
            d += data[i].abs() as u64 * w as u64;
            i += 1;
        }
    }

    T::from_u64(d).unwrap()
}

fn manhattan_impl<'py, T>(
    py: Python<'py>,
    indptr: &[i32],
    indices: &[i32],
    data: &[i32],
    weights: Option<&[i32]>,
    output: &mut [T],
    threads: usize,
) -> PyResult<()> 
where
    T: Float + Element + FromPrimitive,
{
    let n = (output.len() as f32 * 2.0).sqrt().ceil() as usize; 

    // cut the result array into several slices so that rayon can safely
    // write the output in parallel (otherwise, a single mutable slice 
    // cannot be shared across multiple threads)
    let subslices = {
        let mut subslice;
        let mut subslices = Vec::new();
        let mut rest = output;
        for i in 1..n {
            (subslice, rest) = rest.split_at_mut(n - i);
            subslices.push(subslice);
        }
        subslices
    };

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .expect("failed to create thread pool");
   
    py.detach(|| {
        pool.install(|| {
            subslices.into_par_iter()
                .enumerate()
                .map(|(px, d_out)| {
                    // let i_next = indptr[px + 1] as usize;
                    for py in px+1..n {
                        d_out[py - px - 1] = manhattan_pair_impl(
                            indptr,
                            indices,
                            data,
                            weights,
                            px,
                            py,
                        );
                    }
                })
                .collect::<()>();
        });
    });

    Ok(())
}

/// Compute pairwise Manhattan distances for a CSR sparse matrix.
#[pyfunction]
#[pyo3(signature = (data, indices, indptr, weights, distances, threads=0))]
pub fn manhattan<'py>(
    py: Python<'py>,
    data: &Bound<'py, PyArray<i32, numpy::Ix1>>,
    indices: &Bound<'py, PyArray<i32, numpy::Ix1>>,
    indptr: &Bound<'py, PyArray<i32, numpy::Ix1>>,
    weights: Option<&Bound<'py, PyArray<i32, numpy::Ix1>>>,
    distances: &Bound<'py, PyAny>,
    threads: usize,
) -> PyResult<()> {

    let indptr_r = indptr.try_readonly()?;
    let indices_r = indices.try_readonly()?;
    let data_r = data.try_readonly()?;
    
    let indptr_s = indptr_r.as_slice()?;
    let indices_s = indices_r.as_slice()?;
    let data_s = data_r.as_slice()?;

    // weights are optional
    let weights_r;
    let weights_s;
    if let Some(w) = weights {
        weights_r = w.try_readonly()?;
        weights_s = Some(weights_r.as_slice()?);
    } else {
        weights_s = None;
    }
    
    let d = distances.as_borrowed();
    if let Ok(d) = <Bound<'py, PyArray::<f64, numpy::Ix1>>>::extract(d) {
        let mut view = d.try_readwrite()?;
        return manhattan_impl(py, indptr_s, indices_s, data_s, weights_s, view.as_slice_mut()?, threads);
    }
    if let Ok(d) = <Bound<'py, PyArray::<f32, numpy::Ix1>>>::extract(d) {
        let mut view = d.try_readwrite()?;
        return manhattan_impl(py, indptr_s, indices_s, data_s, weights_s, view.as_slice_mut()?, threads);
    }
    if let Ok(d) = <Bound<'py, PyArray::<f16, numpy::Ix1>>>::extract(d) {
        let mut view = d.try_readwrite()?;
        return manhattan_impl(py, indptr_s, indices_s, data_s, weights_s, view.as_slice_mut()?, threads);
    }

    Err(PyTypeError::new_err("Unsupported dtype in `manhattan`"))
}

/// Compute pairwise Manhattan distances for a CSR sparse matrix.
#[pyfunction]
#[pyo3(signature = (data, indices, indptr, weights, i, j))]
pub fn manhattan_pair<'py>(
    _py: Python<'py>,
    data: &Bound<'py, PyArray<i32, numpy::Ix1>>,
    indices: &Bound<'py, PyArray<i32, numpy::Ix1>>,
    indptr: &Bound<'py, PyArray<i32, numpy::Ix1>>,
    weights: Option<&Bound<'py, PyArray<i32, numpy::Ix1>>>,
    i: usize,
    j: usize,
) -> PyResult<f64> {
    let indptr_r = indptr.try_readonly()?;
    let indices_r = indices.try_readonly()?;
    let data_r = data.try_readonly()?;
    
    let indptr_s = indptr_r.as_slice()?;
    let indices_s = indices_r.as_slice()?;
    let data_s = data_r.as_slice()?;
    
    // weights are optional
    let weights_r;
    let weights_s;
    if let Some(w) = weights {
        weights_r = w.try_readonly()?;
        weights_s = Some(weights_r.as_slice()?);
    } else {
        weights_s = None;
    }

    return Ok(manhattan_pair_impl(indptr_s, indices_s, data_s, weights_s, i, j));
}