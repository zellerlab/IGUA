extern crate kodama;
extern crate pyo3;
extern crate numpy;
extern crate rayon;

use pyo3::prelude::*;

mod clustering;
mod distance;

#[pymodule]
#[pyo3(name = "hca")]
pub fn init<'py>(_py: Python<'py>, m: &Bound<'py, PyModule>) -> PyResult<()> {
    m.add("__package__", "igua")?;
    m.add_function(wrap_pyfunction!(distance::manhattan, m)?)?;
    m.add_function(wrap_pyfunction!(distance::manhattan_pair, m)?)?;
    m.add_function(wrap_pyfunction!(clustering::linkage, m)?)?;
    Ok(())
}
