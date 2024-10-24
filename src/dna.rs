use strum_macros::EnumIter;

#[allow(clippy::upper_case_acronyms)]
#[derive(Debug, EnumIter, Hash, Clone, Eq, PartialEq)]
pub enum DNA {
    A,
    T,
    C,
    G,
    N,
}
