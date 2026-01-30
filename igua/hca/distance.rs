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

fn manhattan_impl<'py, T>(
    py: Python<'py>,
    indptr: &[i32],
    indices: &[i32],
    data: &[i32],
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
                    let i_next = indptr[px + 1] as usize;
                    for py in px+1..n {
                        let j_next = indptr[py + 1] as usize;

                        let mut i = indptr[px] as usize;
                        let mut j = indptr[py] as usize;
                        let mut d = 0;

                        while i < i_next && j < j_next {
                            match indices[i].cmp(&indices[j]) {
                                Ordering::Equal => {
                                    d += data[i].abs_diff(data[j]) as u64;
                                    i += 1;
                                    j += 1;
                                }
                                Ordering::Less => {
                                    d += data[i].abs() as u64;
                                    i += 1;
                                }
                                Ordering::Greater => {
                                    d += data[j].abs() as u64;
                                    j += 1;
                                }
                            }
                        }

                        if i == i_next {
                            while j < j_next {
                                d += data[j].abs() as u64;
                                j += 1;
                            }
                        } else {
                            while i < i_next {
                                d += data[i].abs() as u64;
                                i += 1;
                            }
                        }

                        d_out[py - px - 1] = T::from_u64(d).unwrap(); // FIXME
                    }
                })
                .collect::<()>();
        });
    });

    Ok(())
}

/// Compute pairwise Manhattan distances for a CSR sparse matrix.
#[pyfunction]
#[pyo3(signature = (data, indices, indptr, distances, threads=0))]
pub fn manhattan<'py>(
    py: Python<'py>,
    data: &Bound<'py, PyArray<i32, numpy::Ix1>>,
    indices: &Bound<'py, PyArray<i32, numpy::Ix1>>,
    indptr: &Bound<'py, PyArray<i32, numpy::Ix1>>,
    distances: &Bound<'py, PyAny>,
    threads: usize,
) -> PyResult<()> {

    let indptr_r = indptr.try_readonly()?;
    let indices_r = indices.try_readonly()?;
    let data_r = data.try_readonly()?;
    
    let indptr_s = indptr_r.as_slice()?;
    let indices_s = indices_r.as_slice()?;
    let data_s = data_r.as_slice()?;
    
    let d = distances.as_borrowed();
    if let Ok(d) = <Bound<'py, PyArray::<f64, numpy::Ix1>>>::extract(d) {
        let mut view = d.try_readwrite()?;
        return manhattan_impl(py, indptr_s, indices_s, data_s, view.as_slice_mut()?, threads);
    }
    if let Ok(d) = <Bound<'py, PyArray::<f32, numpy::Ix1>>>::extract(d) {
        let mut view = d.try_readwrite()?;
        return manhattan_impl(py, indptr_s, indices_s, data_s, view.as_slice_mut()?, threads);
    }
    if let Ok(d) = <Bound<'py, PyArray::<f16, numpy::Ix1>>>::extract(d) {
        let mut view = d.try_readwrite()?;
        return manhattan_impl(py, indptr_s, indices_s, data_s, view.as_slice_mut()?, threads);
    }

    Err(PyTypeError::new_err("Unsupported dtype in `manhattan`"))
}

fn manhattan_pair_impl(
    indptr: &[i32],
    indices: &[i32],
    data: &[i32],
    px: usize,
    py: usize,
) -> PyResult<f64> {
    let i_next = indptr[px + 1] as usize;
    let j_next = indptr[py + 1] as usize;

    let mut i = indptr[px] as usize;
    let mut j = indptr[py] as usize;
    let mut d = 0;

    while i < i_next && j < j_next {
        match indices[i].cmp(&indices[j]) {
            Ordering::Equal => {
                d += data[i].abs_diff(data[j]) as u64;
                i += 1;
                j += 1;
            }
            Ordering::Less => {
                d += data[i].abs() as u64;
                i += 1;
            }
            Ordering::Greater => {
                d += data[j].abs() as u64;
                j += 1;
            }
        }
    }

    if i == i_next {
        while j < j_next {
            d += data[j].abs() as u64;
            j += 1;
        }
    } else {
        while i < i_next {
            d += data[i].abs() as u64;
            i += 1;
        }
    }

    // println!("{px} {py} {d}");
    Ok(d as f64) // FIXME
}

/// Compute pairwise Manhattan distances for a CSR sparse matrix.
#[pyfunction]
#[pyo3(signature = (data, indices, indptr, i, j))]
pub fn manhattan_pair<'py>(
    _py: Python<'py>,
    data: &Bound<'py, PyArray<i32, numpy::Ix1>>,
    indices: &Bound<'py, PyArray<i32, numpy::Ix1>>,
    indptr: &Bound<'py, PyArray<i32, numpy::Ix1>>,
    i: usize,
    j: usize,
) -> PyResult<f64> {

    let indptr_r = indptr.try_readonly()?;
    let indices_r = indices.try_readonly()?;
    let data_r = data.try_readonly()?;
    
    let indptr_s = indptr_r.as_slice()?;
    let indices_s = indices_r.as_slice()?;
    let data_s = data_r.as_slice()?;

    return manhattan_pair_impl(indptr_s, indices_s, data_s, i, j);
}