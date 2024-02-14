//! Python bindings to the rustyms library.

use std::fmt::Debug;

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

use rustyms;
use rustyms::Chemical;

/// Element.
///
/// A chemical element, with its isotopes and their properties.
///
/// Parameters
/// ----------
/// symbol : str
///
#[pyclass]
#[derive(Debug)]
pub struct Element(rustyms::Element);

#[pymethods]
impl Element {
    #[new]
    fn new(symbol: &str) -> PyResult<Self> {
        match rustyms::Element::try_from(symbol) {
            Ok(element) => Ok(Element(element)),
            Err(_) => Err(PyValueError::new_err("Invalid element symbol.")),
        }
    }

    fn __repr__(&self) -> String {
        format!("Element('{}')", self.0.to_string())
    }

    fn __str__(&self) -> String {
        self.0.to_string()
    }

    /// Get all available isotopes (N, mass, abundance).
    ///
    /// Returns
    /// -------
    /// list[tuple[int, float, float]]
    ///
    fn isotopes(&self) -> Vec<(u16, f64, f64)> {
        self.0
            .isotopes()
            .iter()
            .map(|i| (i.0, i.1.value, i.2))
            .collect()
    }

    /// The mass of the specified isotope of this element (if that isotope exists).
    ///
    /// Parameters
    /// ----------
    /// isotope : int
    ///    The isotope number.
    ///
    /// Returns
    /// -------
    /// float | None
    ///
    fn mass(&self, isotope: u16) -> Option<f64> {
        match self.0.mass(isotope) {
            Some(mass) => Some(mass.value),
            None => None,
        }
    }

    /// The average weight of the specified isotope of this element (if that isotope exists).
    ///
    /// Parameters
    /// ----------
    /// isotope : int
    ///     The isotope number.
    ///
    /// Returns
    /// -------
    /// float
    ///
    fn average_weight(&self, isotope: u16) -> Option<f64> {
        match self.0.average_weight(isotope) {
            Some(mass) => Some(mass.value),
            None => None,
        }
    }

    // TODO: Should first be fixed upstream before exposing.
    // /// Gives the most abundant mass based on the number of this isotope
    // ///
    // /// Parameters
    // /// ----------
    // /// n : int
    // ///   The number of atoms of this element. // //[TODO] CORRECT?
    // /// isotope : int
    // ///   The isotope number.
    // ///
    // /// Returns
    // /// -------
    // /// float | None
    // ///
    // fn most_abundant_mass(&self, n: i16, isotope: u16) -> Option<f64> {
    //     match self.0.most_abundant_mass(n, isotope) {
    //         Some(mass) => Some(mass.value),
    //         None => None,
    //     }
    // }
}

impl std::fmt::Display for Element {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{:?}", self)
    }
}

/// Molecular formula.
///
/// A molecular formula: a selection of elements of specified isotopes together forming a structure.
///
#[pyclass]
#[derive(Debug)]
pub struct MolecularFormula(rustyms::MolecularFormula);

#[pymethods]
impl MolecularFormula {
    // #[new]
    // fn new(formula: &str) -> PyResult<Self> {
    //     todo!()
    // }

    fn __repr__(&self) -> String {
        format!("MolecularFormula('{}')", self.0.to_string())
    }

    fn __str__(&self) -> String {
        self.0.to_string()
    }

    /// Add the given element to this formula (while keeping it ordered and simplified)
    ///
    /// Parameters
    /// ----------
    /// element : Element
    ///     The element to add.
    /// isotope : int
    ///     The isotope number of the element to add.
    /// n : int
    ///     The number of atoms of this element to add.
    ///
    fn add(&mut self, element: &Element, isotope: u16, n: i16) {
        self.0.add((element.0.clone(), isotope, n));
    }

    /// Get the elements making this formula.
    fn elements(&self) -> Vec<(Element, u16, i16)> {
        self.0
            .elements()
            .iter()
            .map(|(e, i, n)| (Element(e.clone()), *i, *n))
            .collect()
    }

    // TODO: Get this working (issues with slice type and pyo3)
    // /// Create a new molecular formula with the given global isotope modifications.
    // fn with_global_isotope_modifications(&self, substitutions: Vec<(Element, u16)>) -> Self {
    //     let substitutions = substitutions
    //         .iter()
    //         .map(|(e, i)| (e.0.clone(), *i))
    //         .collect::<Vec<_>>()
    //         .as_slice();
    //     MolecularFormula(self.0.with_global_isotope_modifications(&substitutions))
    // }

    /// Get the number of electrons (the only charged species, any ionic species is saved as that element +/- the correct number of electrons). The inverse of that number is given as the charge.
    ///
    /// Returns
    /// -------
    /// int
    ///
    fn charge(&self) -> i16 {
        self.0.charge()
    }

    /// The mass of the molecular formula of this element, if all element species (isotopes) exists
    ///
    /// Returns
    /// -------
    /// float | None
    ///
    fn monoisotopic_mass(&self) -> Option<f64> {
        match self.0.monoisotopic_mass() {
            Some(mass) => Some(mass.value),
            None => None,
        }
    }

    /// The average weight of the molecular formula of this element, if all element species (isotopes) exists.
    ///
    /// Returns
    /// -------
    /// float | None
    ///
    fn average_weight(&self) -> Option<f64> {
        match self.0.average_weight() {
            Some(mass) => Some(mass.value),
            None => None,
        }
    }

    /// The most abundant mass, meaning the isotope that will have the highest intensity.
    ///
    /// Returns
    /// -------
    /// float | None
    ///
    fn most_abundant_mass(&self) -> Option<f64> {
        match self.0.most_abundant_mass() {
            Some(mass) => Some(mass.value),
            None => None,
        }
    }

    /// Get the mass in the given mode.
    ///
    /// Parameters
    /// ----------
    /// mode : str
    ///    The mode to get the mass in. One of: ``monoisotopic``, ``average``, ``most_abundant``
    ///
    /// Returns
    /// -------
    /// float | None
    ///
    /// Raises
    /// ------
    /// ValueError
    ///   If the mode is not one of the valid modes.
    ///
    fn mass(&self, mode: &str) -> PyResult<Option<f64>> {
        match mode {
            "monoisotopic" => Ok(self.monoisotopic_mass()),
            "average" => Ok(self.average_weight()),
            "most_abundant" => Ok(self.most_abundant_mass()),
            _ => Err(PyValueError::new_err(
                "Invalid mass mode. Should be one of: monoisotopic, average, most_abundant",
            )),
        }
    }

    /// Create a Hill notation from this collections of elements merged with the pro forma notation for specific isotopes.
    ///
    /// Returns
    /// -------
    /// str
    fn hill_notation(&self) -> String {
        self.0.hill_notation()
    }

    /// Create a Hill notation from this collections of elements merged with the pro forma notation for specific isotopes. Using fancy unicode characters for subscript and superscript numbers.
    ///
    /// Returns
    /// -------
    /// str
    fn hill_notation_fancy(&self) -> String {
        self.0.hill_notation_fancy()
    }

    /// Create a Hill notation from this collections of elements encoded in HTML.
    ///
    /// Returns
    /// -------
    /// str
    fn hill_notation_html(&self) -> String {
        self.0.hill_notation_html()
    }
}

/// Amino acid.
///
/// Parameters
/// ----------
/// name : str
///    The name of the amino acid.
///
#[pyclass]
#[derive(Debug)]
pub struct AminoAcid(rustyms::AminoAcid);

#[pymethods]
impl AminoAcid {
    #[new]
    fn new(name: &str) -> PyResult<Self> {
        match rustyms::AminoAcid::try_from(name) {
            Ok(aa) => Ok(AminoAcid(aa)),
            Err(_) => Err(PyValueError::new_err("Invalid amino acid")),
        }
    }

    fn __str__(&self) -> String {
        self.0.char().to_string()
    }

    fn __repr__(&self) -> String {
        self.to_string()
    }

    /// Molecular formula of the amino acid.
    ///
    /// Returns
    /// -------
    /// MolecularFormula
    ///
    fn formula(&self) -> MolecularFormula {
        MolecularFormula(self.0.formula())
    }

    /// Monoisotopic mass of the amino acid.
    ///
    /// Returns
    /// -------
    /// float | None
    ///
    fn monoisotopic_mass(&self) -> Option<f64> {
        match self.0.formula().monoisotopic_mass() {
            Some(mass) => Some(mass.value),
            None => None,
        }
    }
}

impl std::fmt::Display for AminoAcid {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{:?}", self)
    }
}

/// Amino acid modification.
///
/// Parameters
/// ----------
/// name : str
///   The name of the modification.
///
#[pyclass]
#[derive(Debug, Clone)]
pub struct Modification(rustyms::Modification);

#[pymethods]
impl Modification {
    #[new]
    fn new(name: &str) -> PyResult<Self> {
        match rustyms::Modification::try_from(name, 0..name.len(), &mut vec![]) {
            Ok(modification) => Ok(Modification(modification.defined().unwrap())),
            Err(_) => Err(PyValueError::new_err("Invalid modification")),
        }
    }

    fn __str__(&self) -> String {
        self.0.to_string()
    }

    fn __repr__(&self) -> String {
        format!("Modification('{}')", self.0.to_string())
    }

    /// Molecular formula of the modification.
    ///
    /// Returns
    /// -------
    /// MolecularFormula
    ///
    fn formula(&self) -> MolecularFormula {
        MolecularFormula(self.0.formula())
    }

    /// Monoisotopic mass of the modification.
    ///
    /// Returns
    /// -------
    /// float | None
    ///
    fn monoisotopic_mass(&self) -> Option<f64> {
        match self.0.formula().monoisotopic_mass() {
            Some(mass) => Some(mass.value),
            None => None,
        }
    }
}

/// Modification with ambiguous localisation.
///
/// Parameters
/// ----------
/// id : int
///     The id to compare be able to find the other locations where this modifications can be placed.
/// modification : Modification
///     The modification itself.
/// localisation_score : float | None
///     If present the localisation score, meaning the chance/ratio for this modification to show up on this exact spot.
/// group : tuple[str, bool] | None
///     If this is a named group contain the name and track if this is the preferred location or not.
///
#[pyclass]
#[derive(Debug)]
pub struct AmbiguousModification(rustyms::modification::AmbiguousModification);

#[pymethods]
impl AmbiguousModification {
    #[new]
    fn new(
        id: usize,
        modification: Modification,
        localisation_score: Option<f64>,
        group: Option<(String, bool)>,
    ) -> Self {
        AmbiguousModification(rustyms::modification::AmbiguousModification {
            id,
            modification: modification.0,
            localisation_score,
            group,
        })
    }

    fn __repr__(&self) -> String {
        format!(
            "AmbiguousModification(id={}, modification={}, localisation_score={}, group={:?})",
            self.0.id,
            self.0.modification,
            self.0.localisation_score.unwrap_or_default(),
            match &self.0.group {
                Some((s, b)) => format!("({:?}, {:?})", s, b),
                None => String::new(),
            }
        )
    }

    /// The id to compare be able to find the other locations where this modifications can be placed.
    ///
    /// Returns
    /// -------
    /// int
    ///
    #[getter]
    fn id(&self) -> usize {
        self.0.id
    }

    /// The modification itself.
    ///
    /// Returns
    /// -------
    /// Modification
    ///
    #[getter]
    fn modification(&self) -> Modification {
        Modification(self.0.modification.clone())
    }

    /// If present the localisation score, meaning the chance/ratio for this modification to show up on this exact spot.
    ///
    /// Returns
    /// -------
    /// float | None
    ///
    #[getter]
    fn localisation_score(&self) -> Option<f64> {
        self.0.localisation_score
    }

    /// If this is a named group contain the name and track if this is the preferred location or not.
    ///
    /// Returns
    /// -------
    /// tuple[str, bool] | None
    ///
    #[getter]
    fn group(&self) -> Option<(String, bool)> {
        match &self.0.group {
            Some((s, b)) => Some((s.to_string(), b.clone())),
            None => None,
        }
    }
}

/// A theoretical fragment of a peptide.
#[pyclass]
#[derive(Debug)]
pub struct Fragment(rustyms::Fragment);

#[pymethods]
impl Fragment {
    fn __repr__(&self) -> String {
        format!(
            "Fragment(formula='{:?}', charge={}, ion='{}', peptide_index={}, neutral_loss='{:?}', label='{}')",
            self.formula(),
            self.charge(),
            self.ion(),
            self.peptide_index(),
            self.neutral_loss(),
            self.label()
        )
    }

    /// The theoretical composition.
    ///
    /// Returns
    /// -------
    /// MolecularFormula
    ///
    #[getter]
    fn formula(&self) -> MolecularFormula {
        MolecularFormula(self.0.formula.clone())
    }

    /// The charge.
    ///
    /// Returns
    /// -------
    /// int
    ///
    #[getter]
    fn charge(&self) -> i16 {
        self.0.charge.value as i16
    }

    /// All possible annotations for this fragment.
    ///
    /// Returns
    /// -------
    /// str
    ///
    #[getter]
    fn ion(&self) -> String {
        self.0.ion.to_string() // TODO: Return as exposed Fragment type
    }

    /// The peptide this fragment comes from, saved as the index into the list of peptides in the overarching crate::ComplexPeptide struct.
    ///
    /// Returns
    /// -------
    /// int
    ///
    #[getter]
    fn peptide_index(&self) -> usize {
        self.0.peptide_index
    }

    /// Any neutral losses applied.
    ///
    /// Returns
    /// -------
    /// str | None
    ///
    #[getter]
    fn neutral_loss(&self) -> Option<String> {
        match &self.0.neutral_loss {
            Some(nl) => Some(nl.to_string()),
            None => None,
        }
    }

    /// Additional description for humans.
    ///
    /// Returns
    /// -------
    /// str
    ///
    #[getter]
    fn label(&self) -> String {
        self.0.label.clone()
    }
}

/// One block in a sequence meaning an amino acid and its accompanying modifications.
#[pyclass]
pub struct SequenceElement(rustyms::SequenceElement);

#[pymethods]
impl SequenceElement {
    fn __repr__(&self) -> String {
        format!("SequenceElement(amino_acid='{}', modifications='{:?}', possible_modifications='{:?}', ambiguous='{:?}')", self.aminoacid(), self.modifications(), self.possible_modifications(), self.ambiguous())
    }

    /// The amino acid.
    ///
    /// Returns
    /// -------
    /// AminoAcid
    ///
    #[getter]
    fn aminoacid(&self) -> AminoAcid {
        AminoAcid(self.0.aminoacid)
    }

    /// All present modifications.
    ///
    /// Returns
    /// -------
    /// list[Modification]
    ///
    #[getter]
    fn modifications(&self) -> Vec<Modification> {
        self.0
            .modifications
            .iter()
            .map(|m| Modification(m.clone()))
            .collect()
    }

    /// All ambiguous modifications (could be placed here or on another position).
    ///
    /// Returns
    /// -------
    /// list[AmbiguousModification]
    ///
    #[getter]
    fn possible_modifications(&self) -> Vec<AmbiguousModification> {
        self.0
            .possible_modifications
            .iter()
            .map(|m| AmbiguousModification(m.clone()))
            .collect()
    }

    /// If this amino acid is part of an ambiguous sequence group `(QA)?` in ProForma
    ///
    /// Returns
    /// -------
    /// int | None
    ///
    #[getter]
    fn ambiguous(&self) -> Option<usize> {
        self.0.ambiguous
    }

    /// Get the molecular formula for this position (unless it is B/Z) with the selected ambiguous modifications, without any global isotype modifications.
    ///
    /// Parameters
    /// ----------
    /// selected_ambiguous : int
    ///
    /// Returns
    /// -------
    /// MolecularFormula
    ///
    fn formula(&self, selected_ambiguous: usize) -> MolecularFormula {
        MolecularFormula(self.0.formula(&[selected_ambiguous]).unwrap())
    }

    /// Get the molecular formula for this position (unless it is B/Z) with the ambiguous modifications placed on the very first placed (and updating this in placed), without any global isotype modifications
    ///
    /// Parameters
    /// ----------
    /// placed : bool
    ///
    /// Returns
    /// -------
    /// MolecularFormula
    ///
    fn formula_greedy(&self, placed: bool) -> MolecularFormula {
        MolecularFormula(self.0.formula_greedy(&mut [placed]).unwrap())
    }

    /// Get the molecular formula for this position (unless it is B/Z) with all ambiguous modifications, without any global isotype modifications
    ///
    /// Returns
    /// -------
    /// MolecularFormula
    ///
    fn formula_all(&self) -> MolecularFormula {
        MolecularFormula(self.0.formula_all().unwrap())
    }
}

/// A peptide with all data as provided by ProForma 2.0.
///
/// Parameters
/// ----------
/// proforma : str
///     The ProForma string.
///
#[pyclass]
#[derive(Clone)]
pub struct LinearPeptide(rustyms::LinearPeptide);

#[pymethods]
impl LinearPeptide {
    /// Create a new peptide from a ProForma string.
    #[new]
    fn new(proforma: &str) -> Self {
        LinearPeptide(
            rustyms::ComplexPeptide::pro_forma(proforma)
                .unwrap()
                .assume_linear(),
        )
    }

    fn __str__(&self) -> String {
        self.0.to_string()
    }

    fn __repr__(&self) -> String {
        format!("LinearPeptide({})", self.0.to_string())
    }

    fn __len__(&self) -> usize {
        self.0.sequence.len()
    }

    /// Labile modifications, which will not be found in the actual spectrum.
    ///
    /// Returns
    /// -------
    /// list[Modification]
    ///
    #[getter]
    fn labile(&self) -> Vec<Modification> {
        self.0
            .labile
            .iter()
            .map(|x| Modification(x.clone()))
            .collect()
    }

    /// N-terminal modification.
    ///
    /// Returns
    /// -------
    /// Modification | None
    ///
    #[getter]
    fn n_term(&self) -> Option<Modification> {
        match &self.0.n_term {
            Some(m) => Some(Modification(m.clone())),
            None => None,
        }
    }

    /// C-terminal modification.
    ///
    /// Returns
    /// -------
    /// Modification | None
    ///
    #[getter]
    fn c_term(&self) -> Option<Modification> {
        match &self.0.c_term {
            Some(m) => Some(Modification(m.clone())),
            None => None,
        }
    }

    /// Sequence of the peptide including modifications.
    ///
    /// Returns
    /// -------
    /// list[SequenceElement]
    ///
    #[getter]
    fn sequence(&self) -> Vec<SequenceElement> {
        self.0
            .sequence
            .iter()
            .map(|x| SequenceElement(x.clone()))
            .collect()
    }

    /// For each ambiguous modification list all possible positions it can be placed on. Indexed by the ambiguous modification id.
    ///
    /// Returns
    /// -------
    /// list[list[int]]
    ///
    #[getter]
    fn ambiguous_modifications(&self) -> Vec<Vec<usize>> {
        self.0.ambiguous_modifications.clone()
    }

    /// Stripped sequence, meaning the sequence without any modifications.
    ///
    /// Returns
    /// -------
    /// str
    ///
    #[getter]
    fn stripped_sequence(&self) -> String {
        self.0.sequence.iter().map(|x| x.aminoacid.char()).collect()
    }

    // TODO: How to get the charge as an integer?
    // /// The adduct ions, if specified.
    // #[getter]
    // fn charge_carriers(&self) -> Vec<> {
    //     todo!()
    // }

    /// Get a copy of the peptide with its sequence reversed.
    ///
    /// Returns
    /// -------
    /// LinearPeptide
    ///
    fn reverse(&self) -> LinearPeptide {
        LinearPeptide(self.0.reverse())
    }

    /// Gives the formula for the whole peptide. With the global isotope modifications applied.
    ///
    /// Returns
    /// -------
    /// MolecularFormula | None
    ///
    fn formula(&self) -> Option<MolecularFormula> {
        match self.0.formula() {
            Some(f) => Some(MolecularFormula(f)),
            None => None,
        }
    }

    /// Generate the theoretical fragments for this peptide, with the given maximal charge of the fragments, and the given model. With the global isotope modifications applied.
    ///
    /// Parameters
    /// ----------
    /// max_charge : int
    ///    The maximal charge of the fragments.
    /// model : str
    ///   The model to use for the fragmentation. One of: ``all``, ``cid_hcd``, ``etcid``, ``etd``, ``ethcd``
    ///
    /// Returns
    /// -------
    /// list[Fragment]
    ///   The theoretical fragments.
    ///
    fn generate_theoretical_fragments(
        &self,
        max_charge: i16,
        model: &str,
        // peptide_index: usize, TODO: Required for linear peptide?
    ) -> PyResult<Vec<Fragment>> {
        let model = match_model(model)?;
        Ok(self
            .0
            .generate_theoretical_fragments(
                rustyms::system::Charge::new::<rustyms::system::e>(max_charge as f64),
                &model,
                0, // TODO: Don't hard code?
            )
            .unwrap()
            .iter()
            .map(|f| Fragment(f.clone()))
            .collect())
    }
}

#[pyclass]
struct RawPeak(rustyms::spectrum::RawPeak);

#[pymethods]
impl RawPeak {
    fn __repr__(&self) -> String {
        format!(
            "RawPeak(charge={}, mz={}, intensity={})",
            self.charge(),
            self.mz(),
            self.intensity()
        )
    }

    /// The charge of the peak.
    ///
    /// Returns
    /// -------
    /// int
    ///
    #[getter]
    fn charge(&self) -> i16 {
        self.0.charge.value as i16
    }

    /// The m/z value of the peak.
    ///
    /// Returns
    /// -------
    /// float
    ///
    #[getter]
    fn mz(&self) -> f64 {
        self.0.mz.value
    }

    /// The intensity of the peak.
    ///
    /// Returns
    /// -------
    /// float
    ///
    #[getter]
    fn intensity(&self) -> f64 {
        self.0.intensity
    }
}

#[pyclass]
/// Represents an annotated peak in a mass spectrometry spectrum.
struct AnnotatedPeak(rustyms::spectrum::AnnotatedPeak);

#[pymethods]
impl AnnotatedPeak {
    fn __repr__(&self) -> String {
        format!(
            "AnnotatedPeak(experimental_mz={}, intensity={}, charge={}, annotation=[{:?}])",
            self.experimental_mz(),
            self.intensity(),
            self.charge(),
            self.annotation(),
        )
    }

    /// The experimental m/z value of the peak.
    ///
    /// Returns
    /// -------
    /// float
    ///
    #[getter]
    fn experimental_mz(&self) -> f64 {
        self.0.experimental_mz.value
    }

    /// The intensity of the peak.
    ///
    /// Returns
    /// -------
    /// float
    ///
    #[getter]
    fn intensity(&self) -> f64 {
        self.0.intensity
    }

    /// The charge of the peak.
    ///
    /// Returns
    /// -------
    /// int
    ///
    #[getter]
    fn charge(&self) -> i16 {
        self.0.charge.value as i16
    }

    /// The annotation of the peak, if present.
    ///
    /// Returns
    /// -------
    /// list[Fragment]
    ///
    #[getter]
    fn annotation(&self) -> Vec<Fragment> {
        self.0
            .annotation
            .iter()
            .map(|x| Fragment(x.clone()))
            .collect()
    }
}

/// A raw spectrum (meaning not annotated yet)
///
/// Parameters
/// ----------
/// title : str
///     The title of the spectrum.
/// num_scans : int
///     The number of scans.
/// rt : float
///     The retention time.
/// precursor_charge : float
///     The found precursor charge.
/// precursor_mass : float
///     The found precursor mass.
/// mz_array : list[float]
///     The m/z values of the peaks.
/// intensity_array : list[float]
///     The intensities of the peaks.
///
/// Returns
/// -------
/// RawSpectrum
///
#[pyclass]
pub struct RawSpectrum(rustyms::RawSpectrum);

#[pymethods]
impl RawSpectrum {
    /// Create a new raw spectrum.
    #[new]
    fn new(
        title: &str,
        num_scans: u64,
        rt: f64,
        precursor_charge: f64,
        precursor_mass: f64,
        mz_array: Vec<f64>,
        intensity_array: Vec<f64>,
    ) -> Self {
        RawSpectrum(rustyms::RawSpectrum {
            title: title.to_string(),
            num_scans,
            rt: rustyms::system::Time::new::<rustyms::system::s>(rt),
            charge: rustyms::system::Charge::new::<rustyms::system::e>(precursor_charge),
            mass: rustyms::system::Mass::new::<rustyms::system::dalton>(precursor_mass),
            spectrum: mz_array
                .into_iter()
                .zip(intensity_array.into_iter())
                .map(|(mz, i)| rustyms::spectrum::RawPeak {
                    charge: rustyms::system::Charge::new::<rustyms::system::e>(1.0),
                    mz: rustyms::system::MassOverCharge::new::<rustyms::system::mz>(mz),
                    intensity: i,
                })
                .collect(),
            intensity: None,
            sequence: None,
            raw_file: None,
            raw_scan_number: None,
            raw_index: None,
            sample: None,
            period: None,
            cycle: None,
            experiment: None,
            controller_type: None,
            controller_number: None,
        })
    }

    // fn __str__(&self) -> String {
    //     todo!()
    // }

    // fn __repr__(&self) -> String {
    //     todo!()
    // }

    /// The title.
    ///
    /// Returns
    /// -------
    /// str
    ///
    #[getter]
    fn title(&self) -> String {
        self.0.title.clone()
    }

    /// The number of scans.
    ///
    /// Returns
    /// -------
    /// int
    ///
    #[getter]
    fn num_scans(&self) -> u64 {
        self.0.num_scans
    }

    /// The retention time.
    ///
    /// Returns
    /// -------
    /// float
    ///
    #[getter]
    fn rt(&self) -> f64 {
        self.0.rt.value
    }

    /// The found precursor charge.
    ///
    /// Returns
    /// -------
    /// float
    #[getter]
    fn charge(&self) -> f64 {
        self.0.charge.value
    }

    /// The found precursor mass.
    ///
    /// Returns
    /// -------
    /// float
    ///
    #[getter]
    fn mass(&self) -> f64 {
        self.0.mass.value
    }

    /// The peaks of which this spectrum consists.
    ///
    /// Returns
    /// -------
    /// list[RawPeak]
    ///
    #[getter]
    fn spectrum(&self) -> Vec<RawPeak> {
        self.0.spectrum.iter().map(|x| RawPeak(x.clone())).collect()
    }

    /// Annotate this spectrum with the given peptide
    ///
    /// Parameters
    /// ----------
    /// peptide : LinearPeptide
    ///   The peptide to annotate the spectrum with.
    /// model : str
    ///  The model to use for the fragmentation. One of: ``all``, ``cid_hcd``, ``etcid``, ``etd``, ``ethcd``
    ///
    /// Returns
    /// -------
    /// AnnotatedSpectrum
    ///  The annotated spectrum.
    ///
    /// Raises
    /// ------
    /// ValueError
    ///  If the model is not one of the valid models.
    ///
    fn annotate(&self, peptide: LinearPeptide, model: &str) -> PyResult<AnnotatedSpectrum> {
        let model = match_model(model)?;
        let fragments = peptide
            .0
            .generate_theoretical_fragments(self.0.charge, &model, 0);
        Ok(AnnotatedSpectrum(self.0.annotate(
            rustyms::ComplexPeptide::from(peptide.0),
            &fragments.unwrap(),
            &model,
            rustyms::MassMode::Monoisotopic,
        )))
    }
}

/// An annotated spectrum.
#[pyclass]
pub struct AnnotatedSpectrum(rustyms::AnnotatedSpectrum);

#[pymethods]
impl AnnotatedSpectrum {
    // fn __str__(&self) -> String {
    //     todo!()
    // }

    // fn __repr__(&self) -> String {
    //     todo!()
    // }

    /// The title.
    ///
    /// Returns
    /// -------
    /// str
    ///
    #[getter]
    fn title(&self) -> String {
        self.0.title.clone()
    }

    /// The number of scans.
    ///
    /// Returns
    /// -------
    /// int
    ///
    #[getter]
    fn num_scans(&self) -> u64 {
        self.0.num_scans
    }

    /// The retention time.
    ///
    /// Returns
    /// -------
    /// float
    ///
    #[getter]
    fn rt(&self) -> f64 {
        self.0.rt.value
    }

    /// The found precursor charge.
    ///
    /// Returns
    /// -------
    /// float
    #[getter]
    fn charge(&self) -> f64 {
        self.0.charge.value
    }

    /// The found precursor mass.
    ///
    /// Returns
    /// -------
    /// float
    ///
    #[getter]
    fn mass(&self) -> f64 {
        self.0.mass.value
    }

    /// The peaks of which this spectrum consists.
    ///
    /// Returns
    /// -------
    /// list[AnnotatedPeak]
    ///
    #[getter]
    fn spectrum(&self) -> Vec<AnnotatedPeak> {
        self.0
            .spectrum
            .iter()
            .map(|x| AnnotatedPeak(x.clone()))
            .collect()
    }
}

/// Helper function to match a model string to a rustyms model.
fn match_model(model: &str) -> PyResult<rustyms::Model> {
    match model {
        "all" => Ok(rustyms::Model::all()),
        "cid_hcd" => Ok(rustyms::Model::cid_hcd()),
        "etcid" => Ok(rustyms::Model::ethcd()),
        "etd" => Ok(rustyms::Model::etd()),
        "ethcd" => Ok(rustyms::Model::ethcd()),
        _ => Err(PyValueError::new_err(
            "Invalid model. Should be one of: 'all', 'cid_hcd', 'etcid', 'etd', 'ethcd'",
        )),
    }
}

/// Python bindings to the rustyms library.
#[pymodule]
#[pyo3(name = "rustyms")]
fn rustyms_py03(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<Element>()?;
    m.add_class::<MolecularFormula>()?;
    m.add_class::<AminoAcid>()?;
    m.add_class::<Modification>()?;
    m.add_class::<AmbiguousModification>()?;
    m.add_class::<Fragment>()?;
    m.add_class::<SequenceElement>()?;
    m.add_class::<LinearPeptide>()?;
    m.add_class::<RawPeak>()?;
    m.add_class::<AnnotatedPeak>()?;
    m.add_class::<RawSpectrum>()?;
    m.add_class::<AnnotatedSpectrum>()?;
    Ok(())
}
