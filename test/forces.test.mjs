import assert from 'node:assert/strict';
import test from 'node:test';

import * as d3 from 'd3-force';
import {
    ensureInitialized,
    forceBundle,
    forceCollide,
    forceLink,
    forceManyBody,
    forceX,
    forceY
} from '../src/js/index.js';

await ensureInitialized();

function cloneNodes(nodes) {
    return nodes.map(node => ({ ...node }));
}

function assertNodesClose(actual, expected, tolerance, message) {
    assert.equal(actual.length, expected.length);
    for (let index = 0; index < actual.length; ++index) {
        for (const property of ['x', 'y', 'vx', 'vy']) {
            const difference = Math.abs(actual[index][property] - expected[index][property]);
            assert.ok(
                difference <= tolerance,
                `${message}: node ${index} ${property} differs by ${difference} (${actual[index][property]} vs ${expected[index][property]})`
            );
        }
    }
}

function disposeSimulation(simulation) {
    for (const name of ['link', 'charge', 'collide', 'x', 'y']) {
        simulation.force(name)?.dispose?.();
    }
    simulation.stop();
}

test('link matches D3 with custom id, distance, strength, and iterations', () => {
    const source = [
        { key: 'a', x: 0, y: 0, vx: 0, vy: 0 },
        { key: 'b', x: 100, y: 30, vx: 0, vy: 0 },
        { key: 'c', x: 180, y: -20, vx: 0, vy: 0 }
    ];
    const linkData = [
        { source: 'a', target: 'b', weight: 1 },
        { source: 'b', target: 'c', weight: 2 }
    ];
    const expected = cloneNodes(source);
    const actual = cloneNodes(source);
    const configure = factory => factory(linkData.map(link => ({ ...link })))
        .id(node => node.key)
        .distance(link => 35 + link.weight * 5)
        .strength(link => 0.25 * link.weight)
        .iterations(2);
    const expectedForce = configure(d3.forceLink);
    const actualForce = configure(forceLink);
    const expectedSimulation = d3.forceSimulation(expected).stop().alphaDecay(0).velocityDecay(0).force('link', expectedForce);
    const actualSimulation = d3.forceSimulation(actual).stop().alphaDecay(0).velocityDecay(0).force('link', actualForce);

    expectedSimulation.tick(3);
    actualSimulation.tick(3);
    assertNodesClose(actual, expected, 5e-5, 'link');
    disposeSimulation(expectedSimulation);
    disposeSimulation(actualSimulation);
});

test('collision matches D3 for variable-radius pair resolution', () => {
    const source = [
        { x: 0, y: 0, vx: 1, vy: 0, radius: 8 },
        { x: 10, y: 2, vx: -1, vy: 0, radius: 13 },
        { x: 100, y: 100, vx: 0, vy: 0, radius: 5 }
    ];
    const expected = cloneNodes(source);
    const actual = cloneNodes(source);
    const expectedForce = d3.forceCollide().radius(node => node.radius).strength(0.4).iterations(2);
    const actualForce = forceCollide().radius(node => node.radius).strength(0.4).iterations(2);
    const expectedSimulation = d3.forceSimulation(expected).stop().alphaDecay(0).velocityDecay(0).force('collide', expectedForce);
    const actualSimulation = d3.forceSimulation(actual).stop().alphaDecay(0).velocityDecay(0).force('collide', actualForce);

    expectedSimulation.tick(3);
    actualSimulation.tick(3);
    assertNodesClose(actual, expected, 5e-5, 'collision');
    disposeSimulation(expectedSimulation);
    disposeSimulation(actualSimulation);
});

function quantile(values, probability) {
    const sorted = [...values].sort((left, right) => left - right);
    const position = (sorted.length - 1) * probability;
    const low = Math.floor(position);
    return sorted[low] + ((sorted[low + 1] ?? sorted[low]) - sorted[low]) * (position - low);
}

function layoutMetrics(nodes, links) {
    const centerX = nodes.reduce((sum, node) => sum + node.x, 0) / nodes.length;
    const centerY = nodes.reduce((sum, node) => sum + node.y, 0) / nodes.length;
    const radial = nodes.map(node => Math.hypot(node.x - centerX, node.y - centerY));
    const edgeLengths = links.map(link => Math.hypot(
        link.source.x - link.target.x,
        link.source.y - link.target.y
    ));
    return {
        radialMedian: quantile(radial, 0.5),
        radialP95: quantile(radial, 0.95),
        edgeMedian: quantile(edgeLengths, 0.5),
        edgeP95: quantile(edgeLengths, 0.95)
    };
}

test('accelerated force stack preserves published layout exactly', () => {
    const count = 384;
    const angle = Math.PI * (3 - Math.sqrt(5));
    const baseNodes = Array.from({ length: count }, (_, index) => {
        const radius = 10 * Math.sqrt(0.5 + index);
        return { id: index, x: radius * Math.cos(index * angle), y: radius * Math.sin(index * angle), vx: 0, vy: 0 };
    });
    const baseLinks = [];
    for (let index = 1; index < count; ++index) {
        baseLinks.push({ source: index, target: Math.floor((index - 1) / 2) });
        if (index > 7) baseLinks.push({ source: index, target: (index * 37) % index });
    }

    function run(accelerated) {
        const nodes = cloneNodes(baseNodes);
        const links = baseLinks.map(link => ({ ...link }));
        const simulation = d3.forceSimulation(nodes).stop();
        if (accelerated) {
            simulation.force('link', forceBundle([
                forceLink(links).id(node => node.id).distance(60),
                forceManyBody().strength(-80).theta(0.9).axes(600, 400, 0.005)
            ]));
        } else {
            simulation
                .force('link', d3.forceLink(links).id(node => node.id).distance(60))
                .force('charge', forceManyBody().strength(-80).theta(0.9))
                .force('collide', d3.forceCollide().radius(13).strength(0.1))
                .force('x', forceX(600).strength(0.005))
                .force('y', forceY(400).strength(0.005));
        }
        simulation.tick(300);
        for (const node of nodes) assert.ok([node.x, node.y, node.vx, node.vy].every(Number.isFinite));
        const result = { nodes, metrics: layoutMetrics(nodes, links) };
        disposeSimulation(simulation);
        return result;
    }

    const expected = run(false);
    const actual = run(true);
    assertNodesClose(actual.nodes, expected.nodes, 0, 'accelerated stack');
    assert.deepEqual(actual.metrics, expected.metrics);
});

test('all forces tolerate an empty simulation', () => {
    const simulation = d3.forceSimulation([]).stop()
        .force('link', forceLink([]))
        .force('charge', forceManyBody())
        .force('collide', forceCollide())
        .force('x', forceX(0))
        .force('y', forceY(0));
    assert.doesNotThrow(() => simulation.tick(2));
    disposeSimulation(simulation);
});
