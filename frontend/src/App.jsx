import React, { useState, useEffect } from 'react';
import axios from 'axios';
import Upload from './components/Upload';
import PipelineProgress from './components/PipelineProgress';
import PipelineConfig from './components/PipelineConfig';
import ResultsExplorer from './components/ResultsExplorer';
import ErrorBoundary from './components/ErrorBoundary';

const API_URL = import.meta.env.VITE_API_URL || 'http://localhost:8000/api';

function App() {
    const [view, setView] = useState('upload'); // upload, pipeline, results
    const [taskId, setTaskId] = useState(null);
    const [fileId, setFileId] = useState(null);
    const [pipelineResults, setPipelineResults] = useState(null);
    const [pipelinePlots, setPipelinePlots] = useState({});
    const [currentStep, setCurrentStep] = useState(null);
    const [completedSteps, setCompletedSteps] = useState([]);
    const [pipelineConfig, setPipelineConfig] = useState({});
    const [isRunning, setIsRunning] = useState(false);

    // Poll for pipeline status
    useEffect(() => {
        let interval;
        if (isRunning && taskId) {
            interval = setInterval(async () => {
                try {
                    const res = await axios.get(`${API_URL}/status/${taskId}`);

                    if (res.data.status === 'SUCCESS') {
                        setPipelineResults(res.data.result);
                        setCompletedSteps(Object.keys(res.data.result?.steps || {}));
                        setCurrentStep(null);
                        setIsRunning(false);
                        setView('results');

                        // Fetch plots
                        try {
                            const plotsRes = await axios.get(`${API_URL}/pipeline/${taskId}/plots`);
                            setPipelinePlots(plotsRes.data.plots || {});
                        } catch (e) {
                            console.log('Could not fetch plots');
                        }
                    } else if (res.data.status === 'FAILURE') {
                        setIsRunning(false);
                        setCurrentStep(null);
                    } else if (res.data.status === 'PROGRESS') {
                        const progress = res.data.progress || {};
                        setCurrentStep(progress.step_name || progress.step);
                    }
                } catch (err) {
                    console.error("Polling error", err);
                }
            }, 3000);
        }
        return () => clearInterval(interval);
    }, [isRunning, taskId]);

    const handleUploadComplete = (data) => {
        setTaskId(data.task_id);
        setFileId(data.file_id);
        setIsRunning(true);
        setCurrentStep('qc');
        setView('pipeline');
    };

    const handleStartFullPipeline = async (config) => {
        if (!fileId) return;

        try {
            const res = await axios.post(`${API_URL}/pipeline/start?file_id=${fileId}`, config);
            setTaskId(res.data.task_id);
            setIsRunning(true);
            setCurrentStep('qc');
            setPipelineConfig(config);
        } catch (err) {
            console.error('Failed to start pipeline:', err);
        }
    };

    const handleQCComplete = (result) => {
        setPipelineResults(result);
        setCompletedSteps(['qc']);
        setPipelinePlots(result.plots || {});
        setIsRunning(false);
        setView('pipeline');
    };

    return (
        <div className="container">
            <header style={{ padding: '2rem 0', textAlign: 'center' }}>
                <h1 style={{
                    fontSize: '2.5rem',
                    marginBottom: '0.5rem',
                    background: 'linear-gradient(to right, #6366f1, #ec4899)',
                    WebkitBackgroundClip: 'text',
                    WebkitTextFillColor: 'transparent'
                }}>
                    🧬 Antigravity scRNA-seq
                </h1>
                <p style={{ color: 'var(--text-secondary)', fontSize: '1.1rem' }}>
                    Multi-Agent Single-Cell Analysis Platform
                </p>
            </header>

            {/* Navigation */}
            <nav style={{
                display: 'flex',
                justifyContent: 'center',
                gap: '1rem',
                marginBottom: '2rem'
            }}>
                {['upload', 'pipeline', 'results'].map((v) => (
                    <button
                        key={v}
                        onClick={() => setView(v)}
                        className={view === v ? 'btn-primary' : ''}
                        style={{
                            padding: '0.75rem 1.5rem',
                            borderRadius: '8px',
                            border: view === v ? 'none' : '1px solid var(--glass-border)',
                            background: view === v ? undefined : 'transparent',
                            color: 'var(--text-primary)',
                            cursor: 'pointer',
                            textTransform: 'capitalize',
                            fontWeight: view === v ? 600 : 400
                        }}
                    >
                        {v === 'upload' && '📤 '}
                        {v === 'pipeline' && '🔄 '}
                        {v === 'results' && '📊 '}
                        {v}
                    </button>
                ))}
            </nav>

            <main>
                <ErrorBoundary>
                    {view === 'upload' && (
                        <Upload
                            onUploadComplete={handleUploadComplete}
                            onQCComplete={handleQCComplete}
                        />
                    )}

                    {view === 'pipeline' && (
                        <>
                            <PipelineProgress
                                currentStep={currentStep}
                                completedSteps={completedSteps}
                                stepResults={pipelineResults?.steps}
                            />

                            {!isRunning && completedSteps.includes('qc') && (
                                <PipelineConfig
                                    onConfigChange={setPipelineConfig}
                                    onStartPipeline={handleStartFullPipeline}
                                />
                            )}

                            {isRunning && (
                                <div className="glass-panel" style={{
                                    padding: '2rem',
                                    textAlign: 'center',
                                    marginTop: '1rem'
                                }}>
                                    <div className="spinner" style={{
                                        fontSize: '2rem',
                                        marginBottom: '1rem',
                                        animation: 'spin 1s linear infinite'
                                    }}>
                                        🧬
                                    </div>
                                    <p>Pipeline running... Processing <strong>{currentStep}</strong></p>
                                    <p style={{ color: 'var(--text-secondary)', fontSize: '0.9rem' }}>
                                        Agents are coordinating analysis steps
                                    </p>
                                </div>
                            )}

                            {pipelineResults && (
                                <div style={{ marginTop: '1rem' }}>
                                    <ResultsExplorer
                                        results={pipelineResults}
                                        plots={pipelinePlots}
                                    />
                                </div>
                            )}
                        </>
                    )}

                    {view === 'results' && (
                        <ResultsExplorer
                            results={pipelineResults}
                            plots={pipelinePlots}
                        />
                    )}
                </ErrorBoundary>
            </main>

            <footer style={{
                marginTop: '4rem',
                textAlign: 'center',
                color: 'var(--text-secondary)',
                fontSize: '0.8rem',
                paddingBottom: '2rem'
            }}>
                <p>© 2025 Antigravity Bio • Powered by Multi-Agent Systems</p>
                <p style={{ marginTop: '0.5rem' }}>
                    🔬 Scanpy • 🧠 7 Coordinated Agents • 📊 Publication-Ready Results
                </p>
            </footer>

            <style>{`
                @keyframes spin {
                    from { transform: rotate(0deg); }
                    to { transform: rotate(360deg); }
                }
            `}</style>
        </div>
    );
}

export default App;
