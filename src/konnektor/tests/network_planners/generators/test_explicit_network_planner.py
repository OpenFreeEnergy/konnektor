# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

import pytest
from konnektor.network_planners.generators.explicit_network import ExplicitNetworkGenerator

def test_explicit_network_planner():


    planner = ExplicitNetworkGenerator()



"""
def test_network_from_names(atom_mapping_basic_test_files):
    ligs = list(atom_mapping_basic_test_files.values())

    requested = [
        ('toluene', '2-naftanol'),
        ('2-methylnaphthalene', '2-naftanol'),
    ]

    network = openfe.setup.ligand_network_planning.generate_network_from_names(
        ligands=ligs,
        names=requested,
        mapper=openfe.LomapAtomMapper(),
    )

    assert len(network.nodes) == len(ligs)
    assert len(network.edges) == 2
    actual_edges = [(e.componentA.name, e.componentB.name)
                    for e in network.edges]
    assert set(requested) == set(actual_edges)


def test_network_from_names_bad_name(atom_mapping_basic_test_files):
    ligs = list(atom_mapping_basic_test_files.values())

    requested = [
        ('hank', '2-naftanol'),
        ('2-methylnaphthalene', '2-naftanol'),
    ]

    with pytest.raises(KeyError, match="Invalid name"):
        _ = openfe.setup.ligand_network_planning.generate_network_from_names(
            ligands=ligs,
            names=requested,
            mapper=openfe.LomapAtomMapper(),
        )


def test_network_from_names_duplicate_name(atom_mapping_basic_test_files):
    ligs = list(atom_mapping_basic_test_files.values())
    ligs = ligs + [ligs[0]]

    requested = [
        ('toluene', '2-naftanol'),
        ('2-methylnaphthalene', '2-naftanol'),
    ]

    with pytest.raises(ValueError, match="Duplicate names"):
        _ = openfe.setup.ligand_network_planning.generate_network_from_names(
            ligands=ligs,
            names=requested,
            mapper=openfe.LomapAtomMapper(),
        )


def test_network_from_indices(atom_mapping_basic_test_files):
    ligs = list(atom_mapping_basic_test_files.values())

    requested = [(0, 1), (2, 3)]

    network = openfe.setup.ligand_network_planning.generate_network_from_indices(
        ligands=ligs,
        indices=requested,
        mapper=openfe.LomapAtomMapper(),
    )

    assert len(network.nodes) == len(ligs)
    assert len(network.edges) == 2

    edges = list(network.edges)
    expected_edges = {(ligs[0], ligs[1]), (ligs[2], ligs[3])}
    actual_edges = {(edges[0].componentA, edges[0].componentB),
                    (edges[1].componentA, edges[1].componentB)}

    assert actual_edges == expected_edges


def test_network_from_indices_indexerror(atom_mapping_basic_test_files):
    ligs = list(atom_mapping_basic_test_files.values())

    requested = [(20, 1), (2, 3)]

    with pytest.raises(IndexError, match="Invalid ligand id"):
        network = openfe.setup.ligand_network_planning.generate_network_from_indices(
            ligands=ligs,
            indices=requested,
            mapper=openfe.LomapAtomMapper(),
        )


@pytest.mark.parametrize('file_fixture, loader', [
    ['orion_network',
     openfe.setup.ligand_network_planning.load_orion_network],
    ['fepplus_network',
     openfe.setup.ligand_network_planning.load_fepplus_network],
])
def test_network_from_external(file_fixture, loader, request,
                               benzene_modifications):

    network_file = request.getfixturevalue(file_fixture)

    network = loader(
        ligands=[l for l in benzene_modifications.values()],
        mapper=openfe.LomapAtomMapper(),
        network_file=network_file,
    )

    expected_edges = {
        (benzene_modifications['benzene'], benzene_modifications['toluene']),
        (benzene_modifications['benzene'], benzene_modifications['phenol']),
        (benzene_modifications['benzene'], benzene_modifications['benzonitrile']),
        (benzene_modifications['benzene'], benzene_modifications['anisole']),
        (benzene_modifications['benzene'], benzene_modifications['styrene']),
        (benzene_modifications['benzene'], benzene_modifications['benzaldehyde']),
    }

    actual_edges = {(e.componentA, e.componentB) for e in list(network.edges)}

    assert len(network.nodes) == 7
    assert len(network.edges) == 6
    assert actual_edges == expected_edges


@pytest.mark.parametrize('file_fixture, loader', [
    ['orion_network',
     openfe.setup.ligand_network_planning.load_orion_network],
    ['fepplus_network',
     openfe.setup.ligand_network_planning.load_fepplus_network],
])
def test_network_from_external_unknown_edge(file_fixture, loader, request,
                                            benzene_modifications):
    network_file = request.getfixturevalue(file_fixture)
    ligs = [l for l in benzene_modifications.values() if l.name != 'phenol']

    with pytest.raises(KeyError, match="Invalid name"):
        network = loader(
            ligands=ligs,
            mapper=openfe.LomapAtomMapper(),
            network_file=network_file,
        )


BAD_ORION_NETWORK = \"""\
# Total number of edges: 6
# ------------------------
benzene >>> toluene
benzene >> phenol
benzene >> benzonitrile
benzene >> anisole
benzene >> styrene
benzene >> benzaldehyde
\"""


def test_bad_orion_network(benzene_modifications, tmpdir):
    with tmpdir.as_cwd():
        with open('bad_orion_net.dat', 'w') as f:
            f.write(BAD_ORION_NETWORK)

        with pytest.raises(KeyError, match="line does not match"):
            network = openfe.setup.ligand_network_planning.load_orion_network(
                ligands=[l for l in benzene_modifications.values()],
                mapper=openfe.LomapAtomMapper(),
                network_file='bad_orion_net.dat',
            )



BAD_EDGES = \"""\
1c91235:9c91235 benzene -> toluene
1c91235:7876633 benzene -> phenol
1c91235:2a51f95 benzene -> benzonitrile
1c91235:efja0bc benzene -> anisole
1c91235:7877722 benzene -> styrene
1c91235:99930cd benzene -> benzaldehyde
\"""


def test_bad_edges_network(benzene_modifications, tmpdir):
    with tmpdir.as_cwd():
        with open('bad_edges.edges', 'w') as f:
            f.write(BAD_EDGES)

        with pytest.raises(KeyError, match="line does not match"):
            network = openfe.setup.ligand_network_planning.load_fepplus_network(
                ligands=[l for l in benzene_modifications.values()],
                mapper=openfe.LomapAtomMapper(),
                network_file='bad_edges.edges',
            )

"""