import React, { useState } from 'react';

const TABS = [
    { id: 'qc', name: 'Quality Control', icon: '📊' },
    { id: 'cluster', name: 'Clustering', icon: '🎯' },
    { id: 'markers', name: 'Markers', icon: '🏷️' },
    { id: 'annotation', name: 'Cell Types', icon: '📝' },
];

export default function ResultsExplorer({ results, plots }) {
    const [activeTab, setActiveTab] = useState('qc');

    if (!results) {
        return (
            <div className="glass-panel" style={{ padding: '2rem', textAlign: 'center' }}>
                <p style={{ color: 'var(--text-secondary)' }}>No results available yet. Run the pipeline first.</p>
            </div>
        );
    }

    const renderQCTab = () => {
        const qcData = results.steps?.qc?.result?.data || results.qc_results?.data;
        const qcPlots = plots?.qc || results.plots || {};

        return (
            <div>
                <h4 style={{ marginTop: 0 }}>QC Metrics Summary</h4>

                {qcData && (
                    <div style={{
                        display: 'grid',
                        gridTemplateColumns: 'repeat(auto-fit, minmax(150px, 1fr))',
                        gap: '1rem',
                        marginBottom: '1.5rem'
                    }}>
                        <MetricCard label="Cells" value={qcData.n_cells?.toLocaleString()} />
                        <MetricCard label="Genes" value={qcData.n_genes?.toLocaleString()} />
                        <MetricCard label="Median Genes/Cell" value={qcData.qc_summary?.median_genes?.toFixed(0)} />
                        <MetricCard label="Median %MT" value={`${qcData.qc_summary?.median_mito_pct?.toFixed(1)}%`} />
                    </div>
                )}

                {qcPlots.qc_violin && (
                    <div style={{ marginBottom: '1.5rem' }}>
                        <h5>QC Violin Plots</h5>
                        <img
                            src={`data:image/png;base64,${qcPlots.qc_violin}`}
                            alt="QC Violin Plots"
                            style={{ maxWidth: '100%', borderRadius: '8px' }}
                        />
                    </div>
                )}

                {qcPlots.qc_scatter && (
                    <div>
                        <h5>Genes vs Counts</h5>
                        <img
                            src={`data:image/png;base64,${qcPlots.qc_scatter}`}
                            alt="QC Scatter Plot"
                            style={{ maxWidth: '100%', borderRadius: '8px' }}
                        />
                    </div>
                )}
            </div>
        );
    };

    const renderClusterTab = () => {
        const clusterData = results.steps?.cluster?.result?.data;
        const umapPlot = plots?.umap || results.steps?.umap?.plot;

        return (
            <div>
                <h4 style={{ marginTop: 0 }}>Clustering Results</h4>

                {clusterData && (
                    <div style={{ marginBottom: '1.5rem' }}>
                        <div style={{
                            display: 'grid',
                            gridTemplateColumns: 'repeat(auto-fit, minmax(150px, 1fr))',
                            gap: '1rem',
                            marginBottom: '1rem'
                        }}>
                            <MetricCard label="Clusters" value={clusterData.n_clusters} />
                            <MetricCard label="Method" value={clusterData.method} />
                            <MetricCard label="Resolution" value={clusterData.resolution} />
                        </div>

                        <h5>Cluster Sizes</h5>
                        <div style={{ display: 'flex', flexWrap: 'wrap', gap: '0.5rem' }}>
                            {Object.entries(clusterData.cluster_sizes || {}).map(([cluster, size]) => (
                                <span
                                    key={cluster}
                                    style={{
                                        padding: '0.25rem 0.75rem',
                                        background: 'var(--bg-secondary)',
                                        borderRadius: '16px',
                                        fontSize: '0.85rem'
                                    }}
                                >
                                    Cluster {cluster}: {size}
                                </span>
                            ))}
                        </div>
                    </div>
                )}

                {umapPlot && (
                    <div>
                        <h5>UMAP Visualization</h5>
                        <img
                            src={`data:image/png;base64,${umapPlot}`}
                            alt="UMAP Clusters"
                            style={{ maxWidth: '100%', borderRadius: '8px' }}
                        />
                    </div>
                )}
            </div>
        );
    };

    const renderMarkersTab = () => {
        const markersData = results.steps?.markers?.result?.data?.markers;

        return (
            <div>
                <h4 style={{ marginTop: 0 }}>Marker Genes</h4>

                {markersData ? (
                    <div>
                        {Object.entries(markersData).map(([cluster, data]) => (
                            <div key={cluster} style={{ marginBottom: '1.5rem' }}>
                                <h5 style={{
                                    background: 'var(--bg-secondary)',
                                    padding: '0.5rem 1rem',
                                    borderRadius: '8px',
                                    display: 'inline-block'
                                }}>
                                    Cluster {cluster}
                                </h5>
                                <div style={{
                                    display: 'flex',
                                    flexWrap: 'wrap',
                                    gap: '0.5rem',
                                    marginTop: '0.5rem'
                                }}>
                                    {data.genes?.slice(0, 10).map((gene, idx) => (
                                        <span
                                            key={gene}
                                            style={{
                                                padding: '0.25rem 0.75rem',
                                                background: `rgba(99, 102, 241, ${0.8 - idx * 0.07})`,
                                                borderRadius: '16px',
                                                fontSize: '0.85rem'
                                            }}
                                        >
                                            {gene}
                                        </span>
                                    ))}
                                </div>
                            </div>
                        ))}
                    </div>
                ) : (
                    <p style={{ color: 'var(--text-secondary)' }}>No marker data available.</p>
                )}
            </div>
        );
    };

    const renderAnnotationTab = () => {
        const annotationData = results.steps?.annotate?.result?.data;

        return (
            <div>
                <h4 style={{ marginTop: 0 }}>Cell Type Annotation</h4>

                {annotationData ? (
                    <div>
                        <p style={{ color: 'var(--text-secondary)', marginBottom: '1rem' }}>
                            Method: {annotationData.method}
                            {annotationData.model && ` (Model: ${annotationData.model})`}
                        </p>

                        {annotationData.cell_types && (
                            <div style={{
                                display: 'grid',
                                gridTemplateColumns: 'repeat(auto-fit, minmax(200px, 1fr))',
                                gap: '1rem'
                            }}>
                                {Object.entries(annotationData.cell_types).map(([cellType, count]) => (
                                    <div
                                        key={cellType}
                                        style={{
                                            padding: '1rem',
                                            background: 'var(--bg-secondary)',
                                            borderRadius: '8px',
                                            display: 'flex',
                                            justifyContent: 'space-between',
                                            alignItems: 'center'
                                        }}
                                    >
                                        <span style={{ fontWeight: 500 }}>{cellType}</span>
                                        <span style={{
                                            color: 'var(--accent-primary)',
                                            fontWeight: 600
                                        }}>
                                            {count}
                                        </span>
                                    </div>
                                ))}
                            </div>
                        )}
                    </div>
                ) : (
                    <p style={{ color: 'var(--text-secondary)' }}>No annotation data available.</p>
                )}
            </div>
        );
    };

    const renderContent = () => {
        switch (activeTab) {
            case 'qc': return renderQCTab();
            case 'cluster': return renderClusterTab();
            case 'markers': return renderMarkersTab();
            case 'annotation': return renderAnnotationTab();
            default: return null;
        }
    };

    return (
        <div className="glass-panel" style={{ padding: '1.5rem' }}>
            <h3 style={{ marginTop: 0, marginBottom: '1rem', display: 'flex', alignItems: 'center', gap: '0.5rem' }}>
                <span>🔍</span> Results Explorer
            </h3>

            {/* Tab Navigation */}
            <div style={{
                display: 'flex',
                gap: '0.5rem',
                marginBottom: '1.5rem',
                borderBottom: '1px solid var(--glass-border)',
                paddingBottom: '0.5rem'
            }}>
                {TABS.map(tab => (
                    <button
                        key={tab.id}
                        onClick={() => setActiveTab(tab.id)}
                        style={{
                            padding: '0.5rem 1rem',
                            background: activeTab === tab.id
                                ? 'linear-gradient(135deg, var(--accent-primary), var(--accent-secondary))'
                                : 'transparent',
                            border: activeTab === tab.id ? 'none' : '1px solid var(--glass-border)',
                            borderRadius: '8px',
                            color: 'var(--text-primary)',
                            cursor: 'pointer',
                            transition: 'all 0.2s',
                            display: 'flex',
                            alignItems: 'center',
                            gap: '0.5rem',
                            fontSize: '0.9rem'
                        }}
                    >
                        <span>{tab.icon}</span>
                        {tab.name}
                    </button>
                ))}
            </div>

            {/* Tab Content */}
            <div style={{ minHeight: '300px' }}>
                {renderContent()}
            </div>
        </div>
    );
}

function MetricCard({ label, value }) {
    return (
        <div style={{
            padding: '1rem',
            background: 'var(--bg-secondary)',
            borderRadius: '8px',
            textAlign: 'center'
        }}>
            <div style={{
                fontSize: '1.5rem',
                fontWeight: 700,
                background: 'linear-gradient(135deg, var(--accent-primary), var(--accent-secondary))',
                WebkitBackgroundClip: 'text',
                WebkitTextFillColor: 'transparent'
            }}>
                {value || '—'}
            </div>
            <div style={{
                fontSize: '0.85rem',
                color: 'var(--text-secondary)',
                marginTop: '0.25rem'
            }}>
                {label}
            </div>
        </div>
    );
}
