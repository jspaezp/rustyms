//! Resolving the attribute sets claimed by a library entry into that entry's attributes.
use context_error::Context;
use indexmap::IndexMap;
use itertools::Itertools;
use mzcv::{Curie, term};

use crate::mzspeclib::{Attribute, Attributes};

/// The attributes of a single library entry or attribute set, keyed by group ID. The [`None`] key
/// holds the ungrouped attributes. Insertion order is the order the attributes appear in the file.
pub(crate) type AttributeGroups = IndexMap<Option<u32>, Vec<(Attribute, Context<'static>)>>;

/// All attribute sets defined in the library header for a single
/// [`EntryType`](crate::mzspeclib::EntryType), keyed by name, in definition order.
pub(crate) type AttributeSets = IndexMap<String, AttributeGroups>;

/// The reserved name of the attribute set that every entry in the library inherits.
pub(crate) const RESERVED_ALL: &str = "all";

/// Resolve the attribute sets claimed by a single library entry into that entry's attributes.
///
/// Per mzSpecLib v1.0 §4.1.11 the precedence order, lowest first, is:
///
/// 1. the reserved `all` set, which every entry inherits without claiming it,
/// 2. each named set the entry claims, in the order the entry claims them (not the order the sets
///    are defined in the header),
/// 3. the entry's own attributes.
///
/// All instances of one term within a single source are inherited and overridden together as an
/// inseparable set, so a later source replaces every earlier instance of the terms it defines
/// instead of adding to them.
///
/// Groups are separate contexts. A set claimed inside a group is inherited into that group
/// (§4.1.12 example 3), and `all` is never inherited into a group, only into the ungrouped section.
/// Groups defined inside an inherited set are appended under fresh group IDs so they cannot collide
/// with the entry's own groups.
pub(crate) fn resolve_attribute_sets(
    own: &AttributeGroups,
    sets: Option<&AttributeSets>,
) -> AttributeGroups {
    let no_sets = AttributeSets::new();
    let sets = sets.unwrap_or(&no_sets);
    let no_attributes = Vec::new();
    let own_ungrouped = own.get(&None).unwrap_or(&no_attributes);

    // The sources for the ungrouped section, lowest precedence first. A claimed set that is not
    // defined in the header is skipped; validating claims is not this function's job.
    let claimed = claimed_set_names(own_ungrouped);
    let sources: Vec<&AttributeGroups> = std::iter::once(RESERVED_ALL)
        .chain(claimed.iter().map(String::as_str))
        .filter_map(|name| sets.get(name))
        .collect();

    let mut ungrouped = IndexMap::new();
    for source in &sources {
        for (accession, instances) in
            by_accession(source.get(&None).unwrap_or(&no_attributes))
        {
            ungrouped.insert(accession, instances);
        }
    }
    for (accession, instances) in by_accession(own_ungrouped) {
        ungrouped.insert(accession, instances);
    }

    let mut resolved = AttributeGroups::new();
    resolved.insert(None, ungrouped.into_values().flatten().collect());

    // Groups defined inside an inherited set come along with the set. They are keyed on the terms
    // they contain so a later source can override a group contributed by an earlier one.
    let mut inherited_groups: IndexMap<Vec<Curie>, Vec<(Attribute, Context<'static>)>> =
        IndexMap::new();
    for source in &sources {
        for instances in source
            .iter()
            .filter_map(|(group, instances)| group.map(|_| instances))
        {
            let signature = instances
                .iter()
                .map(|(attribute, _)| attribute.name.accession)
                .sorted()
                .collect();
            inherited_groups.insert(signature, instances.clone());
        }
    }

    // The entry's own groups, each resolved against the sets claimed inside that group.
    for (group, instances) in own
        .iter()
        .filter(|(group, _)| group.is_some())
    {
        let mut group_attributes = IndexMap::new();
        for name in claimed_set_names(instances) {
            if let Some(set) = sets.get(&name) {
                for (accession, inherited) in
                    by_accession(set.get(&None).unwrap_or(&no_attributes))
                {
                    group_attributes.insert(accession, inherited);
                }
            }
        }
        for (accession, instances) in by_accession(instances) {
            group_attributes.insert(accession, instances);
        }
        resolved.insert(*group, group_attributes.into_values().flatten().collect());
    }

    let first_free_group = own.keys().flatten().copied().max().map_or(0, |g| g + 1);
    for (group, instances) in (first_free_group..).zip(inherited_groups.into_values()) {
        resolved.insert(Some(group), instances);
    }

    resolved
}

/// Flatten resolved attributes into the [`Attributes`] layout, where index 0 is the ungrouped
/// section and every index after that is one group.
pub(crate) fn flatten_attribute_groups(groups: &AttributeGroups) -> Attributes {
    let mut attributes = vec![
        groups
            .get(&None)
            .into_iter()
            .flatten()
            .map(|(attribute, _)| attribute.clone())
            .collect(),
    ];
    attributes.extend(
        groups
            .iter()
            .filter(|(group, _)| group.is_some())
            .map(|(_, instances)| {
                instances.iter().map(|(attribute, _)| attribute.clone()).collect()
            }),
    );
    attributes
}

/// Group attributes by accession, preserving first-appearance order, so that all instances of one
/// term can be inherited and overridden together.
fn by_accession(
    attributes: &[(Attribute, Context<'static>)],
) -> IndexMap<Curie, Vec<(Attribute, Context<'static>)>> {
    let mut by_accession: IndexMap<_, Vec<_>> = IndexMap::new();
    for attribute in attributes {
        by_accession
            .entry(attribute.0.name.accession)
            .or_default()
            .push(attribute.clone());
    }
    by_accession
}

/// The names of the attribute sets claimed in this context, in the order they are claimed.
fn claimed_set_names(attributes: &[(Attribute, Context<'static>)]) -> Vec<String> {
    attributes
        .iter()
        .filter(|(attribute, _)| {
            attribute.name == term!(MS:1003212|library attribute set name)
        })
        .map(|(attribute, _)| attribute.value.to_string())
        .collect()
}

#[cfg(test)]
#[expect(clippy::missing_panics_doc)]
mod tests {
    use super::*;

    /// Build attribute groups from mzSpecLib text lines, as they would be read from a file.
    fn groups(lines: &[&str]) -> AttributeGroups {
        let mut groups = AttributeGroups::new();
        for line in lines {
            let (group_id, attribute, _) = Attribute::parse(line).unwrap();
            groups.entry(group_id).or_default().push((attribute, Context::default().to_owned()));
        }
        groups
    }

    /// The resolved ungrouped attributes, rendered back to `name=value` text.
    fn ungrouped(resolved: &AttributeGroups) -> Vec<String> {
        resolved
            .get(&None)
            .into_iter()
            .flatten()
            .map(|(attribute, _)| attribute.to_string())
            .collect()
    }

    /// The attribute sets from specification v1.0 §4.1.12 example 1.
    fn example_1_sets() -> AttributeSets {
        let mut sets = AttributeSets::new();
        sets.insert(
            RESERVED_ALL.to_string(),
            groups(&[
                "MS:1000465|scan polarity=MS:1000130|positive scan",
                "MS:1000044|dissociation method=MS:1000133|collision-induced dissociation",
                "MS:1003072|spectrum origin type=MS:1003073|observed spectrum",
            ]),
        );
        sets.insert(
            "ETD".to_string(),
            groups(&[
                "MS:1000031|instrument model=MS:1000639|LTQ Orbitrap XL ETD",
                "MS:1000044|dissociation method=MS:1000598|electron transfer dissociation",
            ]),
        );
        sets.insert(
            "HCD".to_string(),
            groups(&[
                "MS:1000031|instrument model=MS:1002416|Orbitrap Fusion",
                "MS:1000044|dissociation method=MS:1000422|beam-type collision-induced dissociation",
                "MS:1000419|collision gas=helium",
                "MS:1000138|normalized collision energy=35",
            ]),
        );
        sets.insert(
            "Decoy".to_string(),
            groups(&[
                "MS:1003072|spectrum origin type=MS:1003195|shuffle-and-reposition decoy spectrum",
            ]),
        );
        sets
    }

    /// §4.1.12 example 1, first spectrum: the entry's own value beats both claimed sets, and the
    /// last-claimed set beats the earlier one and `all`.
    #[test]
    fn entry_attributes_win_over_claimed_sets() {
        let resolved = resolve_attribute_sets(
            &groups(&[
                "MS:1003212|library attribute set name=ETD",
                "MS:1003212|library attribute set name=Decoy",
                "MS:1003072|spectrum origin type=MS:1003194|precursor shift decoy spectrum",
            ]),
            Some(&example_1_sets()),
        );

        let resolved = ungrouped(&resolved);
        assert!(
            resolved.contains(
                &"MS:1000465|scan polarity=MS:1000130|positive scan".to_string()
            ),
            "inherited from all without conflict: {resolved:?}"
        );
        assert!(
            resolved.contains(
                &"MS:1000044|dissociation method=MS:1000598|electron transfer dissociation"
                    .to_string()
            ),
            "ETD overrides the collision-induced dissociation from all: {resolved:?}"
        );
        assert!(
            resolved.contains(
                &"MS:1003072|spectrum origin type=MS:1003194|precursor shift decoy spectrum"
                    .to_string()
            ),
            "the entry's own origin type overrides Decoy, which overrides all: {resolved:?}"
        );
        assert_eq!(
            resolved
                .iter()
                .filter(|a| a.starts_with("MS:1003072|"))
                .count(),
            1,
            "origin type is not repeatable: {resolved:?}"
        );
    }

    /// §4.1.12 example 1, second spectrum: claim order decides, not header definition order, and a
    /// term only present in an earlier set survives.
    #[test]
    fn claim_order_decides_precedence() {
        let resolved = ungrouped(&resolve_attribute_sets(
            &groups(&[
                "MS:1003212|library attribute set name=HCD",
                "MS:1003212|library attribute set name=ETD",
            ]),
            Some(&example_1_sets()),
        ));

        for expected in [
            "MS:1000031|instrument model=MS:1000639|LTQ Orbitrap XL ETD",
            "MS:1000044|dissociation method=MS:1000598|electron transfer dissociation",
            "MS:1000419|collision gas=helium",
            "MS:1000138|normalized collision energy=35",
            "MS:1000465|scan polarity=MS:1000130|positive scan",
            "MS:1003072|spectrum origin type=MS:1003073|observed spectrum",
        ] {
            assert!(
                resolved.contains(&expected.to_string()),
                "missing {expected}: {resolved:?}"
            );
        }
    }

    /// The reserved `all` set is inherited by an entry that claims nothing at all.
    #[test]
    fn all_applies_without_being_claimed() {
        let resolved = ungrouped(&resolve_attribute_sets(
            &groups(&["MS:1003061|spectrum name=PEPTIDEK/2"]),
            Some(&example_1_sets()),
        ));

        assert!(
            resolved.contains(
                &"MS:1003072|spectrum origin type=MS:1003073|observed spectrum".to_string()
            ),
            "{resolved:?}"
        );
    }

    /// §4.1.12 example 2: multiple instances of one term are overridden together as one set.
    #[test]
    fn repeated_instances_are_overridden_together() {
        let mut sets = AttributeSets::new();
        sets.insert(
            RESERVED_ALL.to_string(),
            groups(&[
                "MS:1000543|data processing action=MS:1000033|deisotoping",
                "MS:1000543|data processing action=MS:1003242|rank transform",
            ]),
        );

        let resolved = ungrouped(&resolve_attribute_sets(
            &groups(&["MS:1000543|data processing action=MS:1003241|square root transform"]),
            Some(&sets),
        ));

        assert_eq!(
            resolved,
            vec!["MS:1000543|data processing action=MS:1003241|square root transform".to_string()],
            "both inherited instances should be replaced, not appended to"
        );
    }

    /// §4.1.12 example 3: a set claimed inside a group is inherited into that group, and `all` is
    /// never inherited into a group.
    #[test]
    fn a_set_claimed_inside_a_group_is_inherited_into_that_group() {
        let mut sets = AttributeSets::new();
        sets.insert(
            RESERVED_ALL.to_string(),
            groups(&["MS:1003072|spectrum origin type=MS:1003073|observed spectrum"]),
        );
        sets.insert(
            "human_tryptic".to_string(),
            groups(&[
                "MS:1001469|taxonomy: scientific name=Homo sapiens",
                "MS:1001045|cleavage agent name=MS:1001251|Trypsin",
                "MS:1003048|number of enzymatic termini=2",
            ]),
        );

        let resolved = resolve_attribute_sets(
            &groups(&[
                "MS:1003270|proforma peptidoform ion notation=DLGEENFK/2",
                "[2]MS:1000885|protein accession=sp|P02768|ALBU_HUMAN",
                "[2]MS:1003212|library attribute set name=human_tryptic",
            ]),
            Some(&sets),
        );

        let group: Vec<_> = resolved[&Some(2)]
            .iter()
            .map(|(attribute, _)| attribute.to_string())
            .collect();
        assert!(
            group.contains(&"MS:1001045|cleavage agent name=MS:1001251|Trypsin".to_string()),
            "the claimed set is inherited into group 2: {group:?}"
        );
        assert!(
            !group.iter().any(|a| a.starts_with("MS:1003072|")),
            "all must not be inherited into a group: {group:?}"
        );
        assert!(
            ungrouped(&resolved)
                .contains(&"MS:1003072|spectrum origin type=MS:1003073|observed spectrum".to_string()),
            "all is still inherited into the ungrouped section"
        );
    }
}
