use context_error::Context;
use indexmap::IndexMap;
use mzcv::{Curie, term};

use crate::mzspeclib::{Attribute, Attributes};

pub(crate) type AttributeGroups = IndexMap<Option<u32>, Vec<(Attribute, Context<'static>)>>;
pub(crate) type AttributeSets = IndexMap<String, AttributeGroups>;

pub(crate) fn resolve_attribute_sets(
    own: &AttributeGroups,
    sets: Option<&AttributeSets>,
) -> AttributeGroups {
    let empty = Vec::new();
    let own_ungrouped = own.get(&None).unwrap_or(&empty);
    let Some(sets) = sets else {
        return own.clone();
    };

    let claimed = claimed_sets(own_ungrouped);
    let sources: Vec<_> = std::iter::once("all")
        .chain(claimed.iter().map(String::as_str))
        .filter_map(|name| sets.get(name))
        .collect();

    let mut resolved = AttributeGroups::new();
    resolved.insert(
        None,
        merge_by_accession(
            sources
                .iter()
                .map(|source| source.get(&None).unwrap_or(&empty))
                .chain(std::iter::once(own_ungrouped))
                .collect(),
        ),
    );

    for (group, attributes) in own.iter().filter(|(group, _)| group.is_some()) {
        let claimed = claimed_sets(attributes);
        resolved.insert(
            *group,
            merge_by_accession(
                claimed
                    .iter()
                    .filter_map(|name| sets.get(name))
                    .map(|source| source.get(&None).unwrap_or(&empty))
                    .chain(std::iter::once(attributes))
                    .collect(),
            ),
        );
    }

    let first_group = own.keys().flatten().max().map_or(0, |group| group + 1);
    let inherited_groups = sources
        .iter()
        .flat_map(|source| source.iter())
        .filter_map(|(group, attributes)| group.map(|_| attributes));
    for (group, attributes) in (first_group..).zip(inherited_groups) {
        resolved.insert(Some(group), attributes.clone());
    }
    resolved
}

pub(crate) fn flatten_attribute_groups(groups: &AttributeGroups) -> Attributes {
    let mut attributes = vec![
        groups
            .get(&None)
            .into_iter()
            .flatten()
            .map(|(attribute, _)| attribute.clone())
            .collect(),
    ];
    attributes.extend(groups.iter().filter_map(|(group, attributes)| {
        group.map(|_| attributes.iter().map(|(attribute, _)| attribute.clone()).collect())
    }));
    attributes
}

fn merge_by_accession(
    sources: Vec<&Vec<(Attribute, Context<'static>)>>,
) -> Vec<(Attribute, Context<'static>)> {
    let mut merged: IndexMap<Curie, Vec<_>> = IndexMap::new();
    for source in sources {
        let mut current: IndexMap<Curie, Vec<_>> = IndexMap::new();
        for attribute in source {
            current
                .entry(attribute.0.name.accession)
                .or_default()
                .push(attribute.clone());
        }
        merged.extend(current);
    }
    merged.into_values().flatten().collect()
}

fn claimed_sets(attributes: &[(Attribute, Context<'static>)]) -> Vec<String> {
    attributes
        .iter()
        .filter(|(attribute, _)| attribute.name == term!(MS:1003212|library attribute set name))
        .map(|(attribute, _)| attribute.value.to_string())
        .collect()
}
