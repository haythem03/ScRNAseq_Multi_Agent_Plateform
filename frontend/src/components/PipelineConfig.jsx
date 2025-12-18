import React, { useState } from 'react';

const DEFAULT_CONFIG = {
    filter: {
        min_genes: 5,
        max_genes: 5000,
        max_mito_pct: 20,
        min_cells_per_gene: 1
    },
    normalize: { method: 'log_normalize' },
    hvg: { n_top_genes: 2000 },
    pca: { n_comps: 50 },
    cluster: { resolution: 1.0, method: 'leiden' },
    markers: { n_genes: 25 }
};

export default function PipelineConfig({ onConfigChange, onStartPipeline }) {
    const [config, setConfig] = useState(DEFAULT_CONFIG);
    const [showAdvanced, setShowAdvanced] = useState(false);

    const handleChange = (section, field, value) => {
        const newConfig = {
            ...config,
            [section]: {
                ...config[section],
                [field]: value
            }
        };
        setConfig(newConfig);
        onConfigChange?.(newConfig);
    };

    return (
        <div className="glass-panel" style={{ padding: '1.5rem', marginBottom: '2rem' }}>
            <h3 style={{ marginTop: 0, marginBottom: '1rem', display: 'flex', alignItems: 'center', gap: '0.5rem' }}>
                <span>⚙️</span> Pipeline Configuration
            </h3>

            {/* Quick Settings */}
            <div style={{
                display: 'grid',
                gridTemplateColumns: 'repeat(auto-fit, minmax(200px, 1fr))',
                gap: '1rem',
                marginBottom: '1rem'
            }}>
                <ConfigField
                    label="Max Mito %"
                    value={config.filter.max_mito_pct}
                    onChange={(v) => handleChange('filter', 'max_mito_pct', Number(v))}
                    type="number"
                    min={0}
                    max={100}
                />
                <ConfigField
                    label="Min Genes"
                    value={config.filter.min_genes}
                    onChange={(v) => handleChange('filter', 'min_genes', Number(v))}
                    type="number"
                    min={0}
                />
                <ConfigField
                    label="Min Cells/Gene"
                    value={config.filter.min_cells_per_gene}
                    onChange={(v) => handleChange('filter', 'min_cells_per_gene', Number(v))}
                    type="number"
                    min={0}
                />
                <ConfigField
                    label="Cluster Resolution"
                    value={config.cluster.resolution}
                    onChange={(v) => handleChange('cluster', 'resolution', Number(v))}
                    type="number"
                    step={0.1}
                    min={0.1}
                    max={3}
                />
            </div>

            {/* Advanced Toggle */}
            <button
                onClick={() => setShowAdvanced(!showAdvanced)}
                style={{
                    background: 'transparent',
                    border: '1px solid var(--glass-border)',
                    color: 'var(--text-secondary)',
                    padding: '0.5rem 1rem',
                    borderRadius: '8px',
                    cursor: 'pointer',
                    marginBottom: '1rem',
                    width: '100%',
                    display: 'flex',
                    justifyContent: 'space-between',
                    alignItems: 'center'
                }}
            >
                <span>Advanced Settings</span>
                <span>{showAdvanced ? '▲' : '▼'}</span>
            </button>

            {/* Advanced Settings */}
            {showAdvanced && (
                <div style={{
                    background: 'rgba(0,0,0,0.2)',
                    padding: '1rem',
                    borderRadius: '8px',
                    marginBottom: '1rem'
                }}>
                    <h5 style={{ marginTop: 0 }}>Filtering</h5>
                    <div style={{ display: 'grid', gridTemplateColumns: 'repeat(2, 1fr)', gap: '1rem', marginBottom: '1rem' }}>
                        <ConfigField
                            label="Max Genes"
                            value={config.filter.max_genes}
                            onChange={(v) => handleChange('filter', 'max_genes', Number(v))}
                            type="number"
                        />
                    </div>

                    <h5>Dimensionality Reduction</h5>
                    <div style={{ display: 'grid', gridTemplateColumns: 'repeat(2, 1fr)', gap: '1rem', marginBottom: '1rem' }}>
                        <ConfigField
                            label="HVGs"
                            value={config.hvg.n_top_genes}
                            onChange={(v) => handleChange('hvg', 'n_top_genes', Number(v))}
                            type="number"
                        />
                        <ConfigField
                            label="PCA Components"
                            value={config.pca.n_comps}
                            onChange={(v) => handleChange('pca', 'n_comps', Number(v))}
                            type="number"
                        />
                    </div>

                    <h5>Clustering</h5>
                    <div style={{ display: 'grid', gridTemplateColumns: 'repeat(2, 1fr)', gap: '1rem' }}>
                        <ConfigField
                            label="Method"
                            value={config.cluster.method}
                            onChange={(v) => handleChange('cluster', 'method', v)}
                            type="select"
                            options={['leiden', 'louvain']}
                        />
                        <ConfigField
                            label="Marker Genes"
                            value={config.markers.n_genes}
                            onChange={(v) => handleChange('markers', 'n_genes', Number(v))}
                            type="number"
                        />
                    </div>
                </div>
            )}

            {/* Start Button */}
            <button
                className="btn-primary"
                style={{ width: '100%' }}
                onClick={() => onStartPipeline?.(config)}
            >
                🚀 Run Full Pipeline
            </button>
        </div>
    );
}

function ConfigField({ label, value, onChange, type = 'text', ...props }) {
    const inputStyle = {
        width: '100%',
        padding: '0.5rem',
        background: 'var(--bg-primary)',
        border: '1px solid var(--glass-border)',
        borderRadius: '6px',
        color: 'var(--text-primary)',
        fontSize: '0.9rem'
    };

    return (
        <div>
            <label style={{
                display: 'block',
                marginBottom: '0.25rem',
                fontSize: '0.85rem',
                color: 'var(--text-secondary)'
            }}>
                {label}
            </label>
            {type === 'select' ? (
                <select
                    value={value}
                    onChange={(e) => onChange(e.target.value)}
                    style={inputStyle}
                >
                    {props.options?.map(opt => (
                        <option key={opt} value={opt}>{opt}</option>
                    ))}
                </select>
            ) : (
                <input
                    type={type}
                    value={value}
                    onChange={(e) => onChange(e.target.value)}
                    style={inputStyle}
                    {...props}
                />
            )}
        </div>
    );
}
