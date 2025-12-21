import React, { useState, useEffect } from 'react';
import axios from 'axios';
import { collectBase64Plots, prettifyPlotLabel } from '../utils/plotUtils';

const API_URL = import.meta.env.VITE_API_URL || 'http://localhost:8000/api';

const TABS = [
    { id: 'qc', name: 'Quality Control', icon: '📊' },
    { id: 'cluster', name: 'Clustering', icon: '🎯' },
    { id: 'markers', name: 'Markers', icon: '🏷️' },
    { id: 'annotation', name: 'Cell Types', icon: '📝' },
    { id: 'files', name: 'Data Files', icon: '📁' },
];

export default function ResultsExplorer({ results, plots }) {
    const [activeTab, setActiveTab] = useState('qc');
    const [checkpoints, setCheckpoints] = useState([]);
    const [visualizationFiles, setVisualizationFiles] = useState([]);
    const [selectedFile, setSelectedFile] = useState(null);
    const [fileViz, setFileViz] = useState(null);
    const [vizLoading, setVizLoading] = useState(false);

    // Fetch available checkpoint files and visualization images
    useEffect(() => {
        axios.get(`${API_URL}/checkpoints`)
            .then(res => setCheckpoints(res.data.checkpoints || []))
            .catch(() => setCheckpoints([]));
        axios.get(`${API_URL}/visualizations`)
            .then(res => setVisualizationFiles(res.data.files || []))
            .catch(() => setVisualizationFiles([]));
    }, [results]);

    if (!results) {
        return (
            <div className="glass-panel" style={{ padding: '2rem', textAlign: 'center' }}>
                <p style={{ color: 'var(--text-secondary)' }}>No results available yet. Run the pipeline first.</p>
            </div>
        );
    }

    const qcPlotsFromResults = collectBase64Plots(results, 'qc');
    const qcPlotsFromProps = collectBase64Plots({ plots }, 'qc');
    const qcPlotEntries = Object.entries({ ...qcPlotsFromResults, ...qcPlotsFromProps });

    // Collect all plots from all steps for easy access
    const allPlots = { ...plots };
    if (results.steps) {
        Object.entries(results.steps).forEach(([stepName, stepData]) => {
            if (stepData.plot) allPlots[stepName] = stepData.plot;
            if (stepData.plots) allPlots[stepName] = stepData.plots;
        });
    }
    if (results.plots) {
        Object.entries(results.plots).forEach(([key, value]) => {
            allPlots[key] = value;
        });
    }

    const renderQCTab = () => {
        const qcData = results.steps?.qc?.result?.data || results.qc_results?.data;

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

                {qcPlotEntries.length > 0 && (
                    <div style={{ marginBottom: '1.5rem' }}>
                        <h5 style={{ marginBottom: '0.75rem' }}>QC Visualizations</h5>
                        <div style={{
                            display: 'grid',
                            gridTemplateColumns: 'repeat(auto-fit, minmax(220px, 1fr))',
                            gap: '1rem'
                        }}>
                            {qcPlotEntries.map(([key, value]) => (
                                <div key={key} style={{
                                    background: 'var(--bg-secondary)',
                                    borderRadius: '10px',
                                    padding: '0.75rem',
                                    boxShadow: '0 4px 10px rgba(0,0,0,0.2)'
                                }}>
                                    <p style={{
                                        margin: '0 0 0.5rem',
                                        fontSize: '0.9rem',
                                        color: 'var(--text-secondary)'
                                    }}>
                                        {prettifyPlotLabel(key)}
                                    </p>
                                    <img
                                        src={`data:image/png;base64,${value}`}
                                        alt={prettifyPlotLabel(key)}
                                        style={{ maxWidth: '100%', borderRadius: '8px' }}
                                    />
                                </div>
                            ))}
                        </div>
                    </div>
                )}

                {!qcData && qcPlotEntries.length === 0 && (
                    <p style={{ color: 'var(--text-secondary)' }}>No QC data available yet. Run the pipeline first.</p>
                )}
            </div>
        );
    };

    const renderClusterTab = () => {
        const clusterData = results.steps?.cluster?.result?.data;
        // Get UMAP plot from multiple possible locations
        const umapPlot = allPlots.umap || results.steps?.umap?.plot || plots?.umap;

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

                {umapPlot ? (
                    <div style={{ marginTop: '1.5rem' }}>
                        <h5>UMAP Visualization</h5>
                        <div style={{
                            background: 'var(--bg-secondary)',
                            borderRadius: '10px',
                            padding: '1rem',
                            boxShadow: '0 4px 10px rgba(0,0,0,0.2)'
                        }}>
                            <img
                                src={`data:image/png;base64,${umapPlot}`}
                                alt="UMAP Clusters"
                                style={{ maxWidth: '100%', borderRadius: '8px' }}
                            />
                        </div>
                    </div>
                ) : (
                    <p style={{ color: 'var(--text-secondary)', marginTop: '1rem' }}>
                        {clusterData ? 'UMAP plot not available yet.' : 'No clustering data available. Run the full pipeline first.'}
                    </p>
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
            case 'files': return renderFilesTab();
            default: return null;
        }
    };

    const handleVisualize = async (filename) => {
        setSelectedFile(filename);
        setVizLoading(true);
        try {
            const res = await axios.get(`${API_URL}/checkpoints/${filename}/visualize`);
            setFileViz(res.data);
        } catch (e) {
            console.error('Failed to visualize file:', e);
            setFileViz({ error: e.message });
        }
        setVizLoading(false);
    };

    const renderFilesTab = () => {

        return (
            <div>
                <h4 style={{ marginTop: 0 }}>Generated Data Files</h4>
                <p style={{ color: 'var(--text-secondary)', marginBottom: '1rem' }}>
                    Checkpoint files saved during the analysis pipeline. Click "Visualize" to see plots generated from the h5ad file.
                </p>

                {checkpoints.length > 0 ? (
                    <div style={{
                        display: 'grid',
                        gridTemplateColumns: 'repeat(auto-fit, minmax(320px, 1fr))',
                        gap: '1rem',
                        marginBottom: '1.5rem'
                    }}>
                        {checkpoints.map((file) => (
                            <div
                                key={file.name}
                                style={{
                                    padding: '1rem',
                                    background: selectedFile === file.name ? 'rgba(99, 102, 241, 0.15)' : 'var(--bg-secondary)',
                                    borderRadius: '10px',
                                    border: selectedFile === file.name ? '1px solid var(--accent-primary)' : '1px solid transparent',
                                    boxShadow: '0 2px 8px rgba(0,0,0,0.15)'
                                }}
                            >
                                <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: '0.5rem' }}>
                                    <div>
                                        <div style={{ fontWeight: 500 }}>📄 {file.name}</div>
                                        <div style={{ fontSize: '0.8rem', color: 'var(--text-secondary)' }}>
                                            {file.size_mb} MB • {file.modified}
                                        </div>
                                    </div>
                                </div>
                                <div style={{ display: 'flex', gap: '0.5rem' }}>
                                    <button
                                        onClick={() => handleVisualize(file.name)}
                                        style={{
                                            padding: '0.4rem 0.8rem',
                                            background: 'linear-gradient(135deg, var(--accent-primary), var(--accent-secondary))',
                                            border: 'none',
                                            borderRadius: '6px',
                                            color: 'white',
                                            fontSize: '0.8rem',
                                            fontWeight: 500,
                                            cursor: 'pointer'
                                        }}
                                    >
                                        🔬 Visualize
                                    </button>
                                    <a
                                        href={`${API_URL}/checkpoints/${file.name}`}
                                        download={file.name}
                                        style={{
                                            padding: '0.4rem 0.8rem',
                                            background: 'var(--bg-primary)',
                                            border: '1px solid var(--glass-border)',
                                            borderRadius: '6px',
                                            color: 'var(--text-primary)',
                                            textDecoration: 'none',
                                            fontSize: '0.8rem',
                                            fontWeight: 500
                                        }}
                                    >
                                        ⬇️ Download
                                    </a>
                                </div>
                            </div>
                        ))}
                    </div>
                ) : (
                    <p style={{ color: 'var(--text-secondary)' }}>
                        No checkpoint files available yet. Run the pipeline to generate data files.
                    </p>
                )}

                {/* PNG Visualization Files */}
                {visualizationFiles.length > 0 && (
                    <div style={{ marginBottom: '2rem' }}>
                        <h5 style={{ marginBottom: '0.75rem' }}>📊 Generated Plots</h5>
                        <div style={{
                            display: 'grid',
                            gridTemplateColumns: 'repeat(auto-fit, minmax(280px, 1fr))',
                            gap: '1rem'
                        }}>
                            {visualizationFiles.map((file) => (
                                <div key={file.name} style={{
                                    background: 'var(--bg-secondary)',
                                    borderRadius: '10px',
                                    padding: '0.75rem',
                                    boxShadow: '0 2px 8px rgba(0,0,0,0.15)'
                                }}>
                                    <p style={{ margin: '0 0 0.5rem', fontSize: '0.85rem', color: 'var(--text-secondary)' }}>
                                        {file.name.replace('.png', '').replace(/_/g, ' ')}
                                    </p>
                                    <img
                                        src={`${API_URL}/visualizations/${file.name}`}
                                        alt={file.name}
                                        style={{ maxWidth: '100%', borderRadius: '8px' }}
                                    />
                                </div>
                            ))}
                        </div>
                    </div>
                )}

                {vizLoading && (
                    <div style={{ textAlign: 'center', padding: '2rem' }}>
                        <div style={{ fontSize: '2rem', animation: 'spin 1s linear infinite' }}>🔬</div>
                        <p>Generating visualizations from {selectedFile}...</p>
                        <style>{`@keyframes spin { from { transform: rotate(0deg); } to { transform: rotate(360deg); } }`}</style>
                    </div>
                )}

                {fileViz && !vizLoading && (
                    <div style={{ marginTop: '1rem', animation: 'fadeIn 0.3s' }}>
                        <h5 style={{ marginBottom: '0.5rem' }}>
                            📊 Visualization: {selectedFile}
                        </h5>
                        
                        {fileViz.error ? (
                            <div style={{ color: 'var(--error)', padding: '1rem', background: 'rgba(239,68,68,0.1)', borderRadius: '8px' }}>
                                Error: {fileViz.error}
                            </div>
                        ) : (
                            <>
                                {/* File Info */}
                                <div style={{
                                    display: 'grid',
                                    gridTemplateColumns: 'repeat(auto-fit, minmax(100px, 1fr))',
                                    gap: '0.75rem',
                                    marginBottom: '1.5rem'
                                }}>
                                    <MetricCard label="Cells" value={fileViz.info?.n_cells?.toLocaleString()} />
                                    <MetricCard label="Genes" value={fileViz.info?.n_genes?.toLocaleString()} />
                                    {fileViz.info?.n_clusters && <MetricCard label="Clusters" value={fileViz.info.n_clusters} />}
                                </div>

                                {/* Plots */}
                                <div style={{
                                    display: 'grid',
                                    gridTemplateColumns: 'repeat(auto-fit, minmax(300px, 1fr))',
                                    gap: '1rem'
                                }}>
                                    {fileViz.plots?.umap && (
                                        <div style={{ background: 'var(--bg-secondary)', borderRadius: '10px', padding: '1rem' }}>
                                            <h6 style={{ margin: '0 0 0.5rem', color: 'var(--text-secondary)' }}>UMAP Embedding</h6>
                                            <img
                                                src={`data:image/png;base64,${fileViz.plots.umap}`}
                                                alt="UMAP"
                                                style={{ maxWidth: '100%', borderRadius: '8px' }}
                                            />
                                        </div>
                                    )}
                                    {fileViz.plots?.pca && (
                                        <div style={{ background: 'var(--bg-secondary)', borderRadius: '10px', padding: '1rem' }}>
                                            <h6 style={{ margin: '0 0 0.5rem', color: 'var(--text-secondary)' }}>PCA</h6>
                                            <img
                                                src={`data:image/png;base64,${fileViz.plots.pca}`}
                                                alt="PCA"
                                                style={{ maxWidth: '100%', borderRadius: '8px' }}
                                            />
                                        </div>
                                    )}
                                    {fileViz.plots?.qc_violin && (
                                        <div style={{ background: 'var(--bg-secondary)', borderRadius: '10px', padding: '1rem', gridColumn: 'span 2' }}>
                                            <h6 style={{ margin: '0 0 0.5rem', color: 'var(--text-secondary)' }}>QC Metrics</h6>
                                            <img
                                                src={`data:image/png;base64,${fileViz.plots.qc_violin}`}
                                                alt="QC Violin"
                                                style={{ maxWidth: '100%', borderRadius: '8px' }}
                                            />
                                        </div>
                                    )}
                                </div>

                                {Object.keys(fileViz.plots || {}).length === 0 && (
                                    <p style={{ color: 'var(--text-secondary)' }}>
                                        No visualizations available for this checkpoint. It may not contain UMAP/PCA embeddings yet.
                                    </p>
                                )}

                                {/* Data Structure Info */}
                                <div style={{ marginTop: '1.5rem' }}>
                                    <h6 style={{ color: 'var(--text-secondary)', marginBottom: '0.5rem' }}>Data Structure</h6>
                                    <div style={{ fontSize: '0.85rem', color: 'var(--text-secondary)' }}>
                                        <div><strong>obs columns:</strong> {fileViz.info?.obs_columns?.join(', ') || 'None'}</div>
                                        <div><strong>obsm keys:</strong> {fileViz.info?.obsm_keys?.join(', ') || 'None'}</div>
                                        <div><strong>layers:</strong> {fileViz.info?.layers?.join(', ') || 'None'}</div>
                                    </div>
                                </div>
                            </>
                        )}
                        <style>{`@keyframes fadeIn { from { opacity: 0; } to { opacity: 1; } }`}</style>
                    </div>
                )}
            </div>
        );
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
