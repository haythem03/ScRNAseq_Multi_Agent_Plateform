import React from 'react';

const PIPELINE_STEPS = [
    { id: 'qc', name: 'Quality Control', icon: '📊' },
    { id: 'filter', name: 'Cell Filtering', icon: '🔬' },
    { id: 'normalize', name: 'Normalization', icon: '📐' },
    { id: 'hvg', name: 'Variable Genes', icon: '🧬' },
    { id: 'pca', name: 'PCA', icon: '📉' },
    { id: 'neighbors', name: 'Neighbors', icon: '🔗' },
    { id: 'cluster', name: 'Clustering', icon: '🎯' },
    { id: 'umap', name: 'UMAP', icon: '🗺️' },
    { id: 'markers', name: 'Markers', icon: '🏷️' },
    { id: 'annotate', name: 'Annotation', icon: '📝' },
];

const clampPercent = (value) => Math.max(0, Math.min(100, Number(value) || 0));

const getExplicitProgressPercent = (stepResults, stepId) => {
    const stepData = stepResults?.[stepId];
    if (!stepData) return null;
    const candidate = stepData.progress ?? stepData.result?.progress;
    if (typeof candidate === 'number') return clampPercent(candidate);
    if (candidate?.percent != null) return clampPercent(candidate.percent);
    if (candidate?.value != null) return clampPercent(candidate.value);
    return null;
};

const getStepProgressPercent = (stepId, currentStep, completedSteps, stepResults) => {
    const explicit = getExplicitProgressPercent(stepResults, stepId);
    if (explicit !== null) return explicit;
    if (completedSteps?.includes(stepId)) return 100;
    if (currentStep === stepId) return 55;
    return 0;
};

export default function PipelineProgress({ currentStep, completedSteps, stepResults }) {
    const getStepStatus = (stepId) => {
        if (completedSteps?.includes(stepId)) return 'completed';
        if (currentStep === stepId) return 'active';
        return 'pending';
    };

    const getStepValidation = (stepId) => {
        const result = stepResults?.[stepId];
        if (!result) return null;
        return result.validation?.decision;
    };

    return (
        <div className="glass-panel" style={{ padding: '1.5rem', marginBottom: '2rem' }}>
            <h3 style={{ marginTop: 0, marginBottom: '1.5rem', display: 'flex', alignItems: 'center', gap: '0.5rem' }}>
                <span>🔄</span> Pipeline Progress
            </h3>

            <div style={{ display: 'flex', flexWrap: 'wrap', gap: '0.5rem' }}>
                {PIPELINE_STEPS.map((step, index) => {
                    const status = getStepStatus(step.id);
                    const validation = getStepValidation(step.id);
                    const progressPercent = getStepProgressPercent(step.id, currentStep, completedSteps, stepResults);

                    return (
                        <div key={step.id} style={{ display: 'flex', alignItems: 'center' }}>
                            <div style={{ display: 'flex', flexDirection: 'column', gap: '0.4rem' }}>
                                <div
                                    style={{
                                        display: 'flex',
                                        alignItems: 'center',
                                        gap: '0.5rem',
                                        padding: '0.5rem 1rem',
                                        borderRadius: '8px',
                                        background: status === 'completed'
                                            ? 'rgba(16, 185, 129, 0.2)'
                                            : status === 'active'
                                                ? 'linear-gradient(135deg, rgba(99, 102, 241, 0.3), rgba(236, 72, 153, 0.3))'
                                                : 'rgba(255, 255, 255, 0.05)',
                                        border: status === 'active'
                                            ? '1px solid var(--accent-primary)'
                                            : '1px solid transparent',
                                        animation: status === 'active' ? 'pulse 2s infinite' : 'none',
                                    }}
                                >
                                    <span style={{ fontSize: '1.2rem' }}>
                                        {status === 'completed' ? '✅' : step.icon}
                                    </span>
                                    <span style={{
                                        fontSize: '0.85rem',
                                        color: status === 'pending' ? 'var(--text-secondary)' : 'var(--text-primary)'
                                    }}>
                                        {step.name}
                                    </span>
                                    {validation === 'STOP' && (
                                        <span title="Validation failed" style={{ color: 'var(--error)' }}>⚠️</span>
                                    )}
                                    {validation === 'PROCEED_WITH_CAUTION' && (
                                        <span title="Proceed with caution" style={{ color: '#fbbf24' }}>⚡</span>
                                    )}
                                </div>
                                <div style={{ width: '150px' }}>
                                    <div style={{
                                        width: '100%',
                                        height: '6px',
                                        background: 'rgba(255, 255, 255, 0.08)',
                                        borderRadius: '999px',
                                        overflow: 'hidden'
                                    }}>
                                        <div style={{
                                            width: `${progressPercent}%`,
                                            height: '100%',
                                            background: status === 'completed'
                                                ? 'var(--success)'
                                                : 'linear-gradient(90deg, var(--accent-primary), var(--accent-secondary))',
                                            transition: 'width 0.3s ease'
                                        }} />
                                    </div>
                                    <div style={{
                                        fontSize: '0.65rem',
                                        color: 'var(--text-secondary)',
                                        marginTop: '0.25rem'
                                    }}>
                                        {status === 'completed'
                                            ? 'Complete'
                                            : progressPercent > 0
                                                ? `${progressPercent}%`
                                                : 'Pending'
                                        }
                                    </div>
                                </div>
                            </div>
                            {index < PIPELINE_STEPS.length - 1 && (
                                <div style={{
                                    width: '20px',
                                    height: '2px',
                                    background: status === 'completed'
                                        ? 'var(--success)'
                                        : 'var(--glass-border)',
                                    margin: '0 0.25rem'
                                }} />
                            )}
                        </div>
                    );
                })}
            </div>

            <style>{`
        @keyframes pulse {
          0%, 100% { opacity: 1; }
          50% { opacity: 0.7; }
        }
      `}</style>
        </div>
    );
}
