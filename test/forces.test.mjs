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

test('many-body matches D3 when theta requires direct body evaluation', () => {
    const source = [
        { x: -83, y: 17, vx: 0, vy: 0, strength: -12 },
        { x: -21, y: -54, vx: 0, vy: 0, strength: -27 },
        { x: 13, y: 91, vx: 0, vy: 0, strength: -43 },
        { x: 67, y: -8, vx: 0, vy: 0, strength: -61 },
        { x: 104, y: 73, vx: 0, vy: 0, strength: -79 },
        { x: 149, y: -96, vx: 0, vy: 0, strength: -101 }
    ];
    const expected = cloneNodes(source);
    const actual = cloneNodes(source);
    const configure = factory => factory()
        .strength(node => node.strength)
        .theta(1e-6);
    const expectedSimulation = d3.forceSimulation(expected).stop()
        .alphaDecay(0)
        .velocityDecay(0)
        .force('charge', configure(d3.forceManyBody));
    const actualSimulation = d3.forceSimulation(actual).stop()
        .alphaDecay(0)
        .velocityDecay(0)
        .force('charge', configure(forceManyBody));

    expectedSimulation.tick(1);
    actualSimulation.tick(1);
    assertNodesClose(actual, expected, 2e-5, 'direct many-body');
    disposeSimulation(expectedSimulation);
    disposeSimulation(actualSimulation);
});

test('Barnes-Hut charge is translation invariant', () => {
    const source = Array.from({ length: 48 }, (_, index) => ({
        x: ((index * 47) % 211) - 100,
        y: ((index * 83) % 197) - 90,
        vx: 0,
        vy: 0
    }));

    function run(offsetX, offsetY) {
        const nodes = source.map(node => ({
            ...node,
            x: node.x + offsetX,
            y: node.y + offsetY
        }));
        const simulation = d3.forceSimulation(nodes).stop()
            .alphaDecay(0)
            .velocityDecay(0)
            .force('charge', forceManyBody().strength(-80).theta(0.9));
        simulation.tick(1);
        disposeSimulation(simulation);
        return nodes;
    }

    const origin = run(0, 0);
    const translated = run(700, -350);
    for (let index = 0; index < origin.length; ++index) {
        assert.ok(Math.abs(origin[index].vx - translated[index].vx) < 2e-5, `translated node ${index} vx`);
        assert.ok(Math.abs(origin[index].vy - translated[index].vy) < 2e-5, `translated node ${index} vy`);
        assert.ok(Math.abs(origin[index].x - (translated[index].x - 700)) < 2e-5, `translated node ${index} x`);
        assert.ok(Math.abs(origin[index].y - (translated[index].y + 350)) < 2e-5, `translated node ${index} y`);
    }
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

test('adaptive cadence reduces work as alpha cools', () => {
    const scheduledAlphas = [];
    const regularAlphas = [];
    const scheduled = {
        _prepare() {},
        _apply(alpha) { scheduledAlphas.push(alpha); },
        _flush() {}
    };
    const regular = {
        _prepare() {},
        _apply(alpha) { regularAlphas.push(alpha); },
        _flush() {}
    };
    const bundle = forceBundle([scheduled, regular]).adaptiveCadence(scheduled, {
        targetImpulse: 0.4,
        maxInterval: 8
    });
    const alphas = [1, 0.3, 0.19, 0.18, 0.1, 0.09, 0.08, 0.07, 0.04];

    for (const alpha of alphas) bundle(alpha);

    assert.deepEqual(scheduledAlphas, [1, 0.3, 0.38, 0.4, 0.32]);
    assert.deepEqual(regularAlphas, alphas);
});

test('adaptive cadence can reuse a computed force field between rebuilds', () => {
    const rebuiltAlphas = [];
    const replayedAlphas = [];
    const scheduled = {
        _prepare() {},
        _apply(alpha) { rebuiltAlphas.push(alpha); },
        _replay(alpha) { replayedAlphas.push(alpha); },
        _flush() {}
    };
    const bundle = forceBundle([scheduled]).adaptiveCadence(scheduled, {
        targetImpulse: 0.4,
        maxInterval: 8,
        reuseLast: true
    });
    const alphas = [1, 0.3, 0.19, 0.18, 0.1, 0.09, 0.08, 0.07, 0.04];

    for (const alpha of alphas) bundle(alpha);

    assert.deepEqual(rebuiltAlphas, [1, 0.3, 0.19, 0.1, 0.04]);
    assert.deepEqual(replayedAlphas, [0.18, 0.09, 0.08, 0.07]);
});

test('adaptive cadence rebuilds immediately after a reheat', () => {
    const rebuiltAlphas = [];
    const replayedAlphas = [];
    const scheduled = {
        _prepare() {},
        _apply(alpha) { rebuiltAlphas.push(alpha); },
        _replay(alpha) { replayedAlphas.push(alpha); },
        _flush() {},
        _preparePersistent() { return false; },
        _integrate() {},
        _flushPersistent() {}
    };
    const bundle = forceBundle([scheduled])
        .adaptiveCadence(scheduled, {
            targetImpulse: 0.4,
            maxInterval: 128,
            reuseLast: true
        })
        .substeps(4);

    const cooledAlpha = 0.01 * Math.pow(1 - 0.015, 4);
    bundle(0.01);
    bundle(cooledAlpha);
    bundle(0.2);

    assert.deepEqual(rebuiltAlphas, [0.01, 0.2]);
    assert.deepEqual(replayedAlphas, [cooledAlpha]);
});

test('adaptive cadence rebuilds once per frame while a node is fixed', () => {
    let rebuilds = 0;
    let replays = 0;
    const scheduled = {
        _prepare() {},
        _apply() { ++rebuilds; },
        _replay() { ++replays; },
        _flush() {},
        _preparePersistent() { return true; },
        _hasFixedConstraints() { return true; },
        _integrate() {},
        _flushPersistent() {}
    };
    const bundle = forceBundle([scheduled])
        .adaptiveCadence(scheduled, {
            targetImpulse: 0.4,
            maxInterval: 128,
            reuseLast: true
        })
        .substeps(4);

    bundle(0.01);
    bundle(0.01 * Math.pow(1 - 0.015, 4));

    assert.equal(rebuilds, 2);
    assert.equal(replays, 6);
});

test('fixed-node substeps cap the inferred alpha target', () => {
    let fixed = false;
    const appliedAlphas = [];
    const force = {
        _prepare() {},
        _apply(alpha) { appliedAlphas.push(alpha); },
        _flush() {},
        _preparePersistent() { return true; },
        _hasFixedConstraints() { return fixed; },
        _integrate() {},
        _flushPersistent() {}
    };
    const bundle = forceBundle([force]).substeps(4, {
        fixedAlphaTarget: 0.1
    });
    const outerRetention = Math.pow(1 - 0.015, 4);
    bundle(0.01);
    fixed = true;
    bundle(0.01 * outerRetention + 0.3 * (1 - outerRetention));

    const expected = 0.1 + (0.01 - 0.1) * Math.pow(1 - 0.015, 4);
    assert.ok(Math.abs(appliedAlphas.at(-1) - expected) < 1e-12);
});

test('batched virtual ticks preserve ordinary D3 positions', () => {
    const source = Array.from({ length: 32 }, (_, index) => ({
        id: index,
        x: Math.cos(index) * (20 + index),
        y: Math.sin(index) * (20 + index),
        vx: 0,
        vy: 0
    }));
    const sourceLinks = Array.from({ length: source.length - 1 }, (_, index) => ({
        source: index,
        target: index + 1
    }));

    function run(substeps) {
        const nodes = cloneNodes(source);
        const links = sourceLinks.map(link => ({ ...link }));
        const alphaDecay = 0.015;
        const bundle = forceBundle([
            forceLink(links).id(node => node.id).distance(40),
            forceManyBody().strength(-30).theta(0.9).axes(100, 80, 0.005)
        ]).substeps(substeps, { alphaDecay, velocityDecay: 0.4 });
        const simulation = d3.forceSimulation(nodes).stop()
            .alphaDecay(1 - Math.pow(1 - alphaDecay, substeps))
            .velocityDecay(0.4)
            .force('link', bundle);
        simulation.tick(12 / substeps);
        simulation.stop();
        bundle.dispose();
        return nodes;
    }

    const ordinary = run(1);
    const batched = run(4);
    for (let index = 0; index < ordinary.length; ++index) {
        assert.ok(
            Math.abs(ordinary[index].x - batched[index].x) < 1e-4,
            `x difference: ${Math.abs(ordinary[index].x - batched[index].x)}`
        );
        assert.ok(
            Math.abs(ordinary[index].y - batched[index].y) < 1e-4,
            `y difference: ${Math.abs(ordinary[index].y - batched[index].y)}`
        );
    }
});

test('batched virtual ticks preserve cooled drag and release dynamics', () => {
    const source = Array.from({ length: 48 }, (_, index) => ({
        id: index,
        x: Math.cos(index * 1.7) * (30 + index),
        y: Math.sin(index * 1.7) * (30 + index),
        vx: 0,
        vy: 0
    }));
    const sourceLinks = Array.from({ length: source.length - 1 }, (_, index) => ({
        source: index + 1,
        target: Math.floor(index / 2)
    }));
    const alphaDecay = 0.015;

    function create(substeps) {
        const nodes = cloneNodes(source);
        const links = sourceLinks.map(link => ({ ...link }));
        const bundle = forceBundle([
            forceLink(links).id(node => node.id).distance(40),
            forceManyBody().strength(-30).theta(0.5).axes(100, 80, 0.005)
        ]).substeps(substeps, { alphaDecay, velocityDecay: 0.4 });
        const simulation = d3.forceSimulation(nodes).stop()
            .alphaDecay(1 - Math.pow(1 - alphaDecay, substeps))
            .velocityDecay(0.4)
            .force('link', bundle);
        return { nodes, bundle, simulation };
    }

    const ordinary = create(1);
    const batched = create(4);
    const positionTolerance = 1.1e-2;
    ordinary.simulation.tick(1000);
    batched.simulation.tick(250);

    const assertPositionsClose = message => {
        for (let index = 0; index < ordinary.nodes.length; ++index) {
            for (const property of ['x', 'y']) {
                const difference = Math.abs(ordinary.nodes[index][property] - batched.nodes[index][property]);
                assert.ok(difference < positionTolerance, `${message}: node ${index} ${property} differs by ${difference}`);
            }
        }
    };
    assertPositionsClose('cooled');

    const ordinaryDragged = ordinary.nodes[0];
    const batchedDragged = batched.nodes[0];
    const ordinaryStart = [ordinaryDragged.x, ordinaryDragged.y];
    const batchedStart = [batchedDragged.x, batchedDragged.y];
    ordinaryDragged.fx = ordinaryDragged.x;
    ordinaryDragged.fy = ordinaryDragged.y;
    batchedDragged.fx = batchedDragged.x;
    batchedDragged.fy = batchedDragged.y;
    ordinary.simulation.alphaTarget(0.3);
    batched.simulation.alphaTarget(0.3);

    for (let frame = 0; frame < 24; ++frame) {
        const dx = 70 + 20 * Math.sin(frame * 0.7);
        const dy = 30 * Math.cos(frame * 0.5);
        ordinaryDragged.x = ordinaryDragged.fx = ordinaryStart[0] + dx;
        ordinaryDragged.y = ordinaryDragged.fy = ordinaryStart[1] + dy;
        batchedDragged.x = batchedDragged.fx = batchedStart[0] + dx;
        batchedDragged.y = batchedDragged.fy = batchedStart[1] + dy;
        ordinary.simulation.tick(4);
        batched.simulation.tick(1);
        assertPositionsClose(`drag frame ${frame}`);
    }

    delete ordinaryDragged.fx;
    delete ordinaryDragged.fy;
    delete batchedDragged.fx;
    delete batchedDragged.fy;
    ordinary.simulation.alphaTarget(0);
    batched.simulation.alphaTarget(0);
    for (let frame = 0; frame < 60; ++frame) {
        ordinary.simulation.tick(4);
        batched.simulation.tick(1);
        assertPositionsClose(`release frame ${frame}`);
    }

    ordinary.simulation.stop();
    batched.simulation.stop();
    ordinary.bundle.dispose();
    batched.bundle.dispose();
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
