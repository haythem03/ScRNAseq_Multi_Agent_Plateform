import React, { useState, useEffect } from 'react';
import axios from 'axios';

const API_URL = import.meta.env.VITE_API_URL || 'http://localhost:8000/api';

export default function Upload({ onUploadComplete, onQCComplete }) {
  const [file, setFile] = useState(null);
  const [status, setStatus] = useState('idle'); // idle, uploading, processing, complete, error
  const [uploadProgress, setUploadProgress] = useState(0);
  const [taskId, setTaskId] = useState(null);
  const [result, setResult] = useState(null);
  const [error, setError] = useState(null);
  const [dragActive, setDragActive] = useState(false);

  const handleFileChange = (e) => {
    const selectedFile = e.target.files[0];
    if (selectedFile) {
      setFile(selectedFile);
      setStatus('idle');
      setError(null);
    }
  };

  const handleDrag = (e) => {
    e.preventDefault();
    e.stopPropagation();
    if (e.type === "dragenter" || e.type === "dragover") {
      setDragActive(true);
    } else if (e.type === "dragleave") {
      setDragActive(false);
    }
  };

  const handleDrop = (e) => {
    e.preventDefault();
    e.stopPropagation();
    setDragActive(false);
    if (e.dataTransfer.files && e.dataTransfer.files[0]) {
      setFile(e.dataTransfer.files[0]);
      setStatus('idle');
      setError(null);
    }
  };

  const handleUpload = async () => {
    if (!file) return;

    setStatus('uploading');
    const formData = new FormData();
    formData.append('file', file);

    try {
      const response = await axios.post(`${API_URL}/upload`, formData, {
        headers: { 'Content-Type': 'multipart/form-data' },
        onUploadProgress: (progressEvent) => {
          const percentCompleted = Math.round((progressEvent.loaded * 100) / progressEvent.total);
          setUploadProgress(percentCompleted);
        }
      });

      setTaskId(response.data.task_id);
      setStatus('processing');
      onUploadComplete?.(response.data);
    } catch (err) {
      console.error(err);
      setError('Upload failed: ' + (err.response?.data?.detail || err.message));
      setStatus('error');
    }
  };

  useEffect(() => {
    let interval;
    if (status === 'processing' && taskId) {
      interval = setInterval(async () => {
        try {
          const res = await axios.get(`${API_URL}/status/${taskId}`);
          if (res.data.status === 'SUCCESS') {
            setResult(res.data.result);
            setStatus('complete');
            onQCComplete?.(res.data.result);
            clearInterval(interval);
          } else if (res.data.status === 'FAILURE') {
            setError('Processing failed.');
            setStatus('error');
            clearInterval(interval);
          }
        } catch (err) {
          console.error("Polling error", err);
        }
      }, 2000);
    }
    return () => clearInterval(interval);
  }, [status, taskId]);

  const formatFileSize = (bytes) => {
    if (bytes < 1024) return bytes + ' B';
    if (bytes < 1024 * 1024) return (bytes / 1024).toFixed(1) + ' KB';
    return (bytes / (1024 * 1024)).toFixed(1) + ' MB';
  };

  return (
    <div className="glass-panel" style={{ padding: '2rem', maxWidth: '700px', margin: '0 auto' }}>
      <h2 style={{ marginTop: 0, marginBottom: '1.5rem', display: 'flex', alignItems: 'center', gap: '0.5rem' }}>
        <span>📤</span> Upload Dataset
      </h2>

      <div
        onDragEnter={handleDrag}
        onDragLeave={handleDrag}
        onDragOver={handleDrag}
        onDrop={handleDrop}
        style={{
          border: `2px dashed ${dragActive ? 'var(--accent-primary)' : 'var(--glass-border)'}`,
          borderRadius: '12px',
          padding: '3rem 2rem',
          textAlign: 'center',
          marginBottom: '1.5rem',
          background: dragActive ? 'rgba(99, 102, 241, 0.1)' : 'rgba(255,255,255,0.02)',
          transition: 'all 0.2s'
        }}
      >
        <input
          type="file"
          onChange={handleFileChange}
          accept=".h5ad,.csv,.h5,.loom,.mtx,.txt"
          style={{ display: 'none' }}
          id="file-upload"
        />

        <div style={{ fontSize: '3rem', marginBottom: '1rem' }}>
          {file ? '📁' : '🧬'}
        </div>

        <label htmlFor="file-upload" className="btn-primary" style={{ display: 'inline-block', cursor: 'pointer' }}>
          {file ? 'Change File' : 'Select File'}
        </label>

        {file && (
          <div style={{ marginTop: '1rem', color: 'var(--text-primary)' }}>
            <strong>{file.name}</strong>
            <span style={{ color: 'var(--text-secondary)', marginLeft: '0.5rem' }}>
              ({formatFileSize(file.size)})
            </span>
          </div>
        )}

        <p style={{ color: 'var(--text-secondary)', marginTop: '1rem', fontSize: '0.9rem' }}>
          Drag & drop or click to select • Supports: .h5ad, .csv, .h5, .loom, .mtx, .txt
        </p>
      </div>

      {status === 'idle' && file && (
        <button className="btn-primary" style={{ width: '100%' }} onClick={handleUpload}>
          🚀 Start Analysis
        </button>
      )}

      {status === 'uploading' && (
        <div>
          <div style={{ display: 'flex', justifyContent: 'space-between', marginBottom: '0.5rem' }}>
            <span>Uploading...</span>
            <span>{uploadProgress}%</span>
          </div>
          <div style={{ height: '8px', width: '100%', background: 'var(--bg-secondary)', borderRadius: '4px' }}>
            <div style={{
              height: '100%',
              width: `${uploadProgress}%`,
              background: 'linear-gradient(90deg, var(--accent-primary), var(--accent-secondary))',
              borderRadius: '4px',
              transition: 'width 0.3s'
            }} />
          </div>
        </div>
      )}

      {status === 'processing' && (
        <div style={{ textAlign: 'center' }}>
          <div style={{
            fontSize: '2rem',
            marginBottom: '1rem',
            animation: 'pulse 1.5s infinite'
          }}>
            🔬
          </div>
          <p style={{ fontWeight: 500 }}>Running Quality Control...</p>
          <p style={{ color: 'var(--text-secondary)', fontSize: '0.9rem' }}>
            Agents are analyzing your data
          </p>
          <style>{`
            @keyframes pulse {
              0%, 100% { transform: scale(1); opacity: 1; }
              50% { transform: scale(1.1); opacity: 0.8; }
            }
          `}</style>
        </div>
      )}

      {status === 'complete' && result && (
        <div style={{ animation: 'fadeIn 0.5s' }}>
          <div style={{
            display: 'flex',
            alignItems: 'center',
            gap: '0.5rem',
            color: 'var(--success)',
            marginBottom: '1rem'
          }}>
            <span style={{ fontSize: '1.5rem' }}>✅</span>
            <h3 style={{ margin: 0 }}>QC Complete</h3>
          </div>

          {result.qc_results?.status === 'error' ? (
            <div style={{
              color: 'var(--error)',
              background: 'rgba(239, 68, 68, 0.1)',
              padding: '1rem',
              borderRadius: '8px'
            }}>
              <strong>Error:</strong> {result.qc_results.message || "Unknown error"}
            </div>
          ) : (
            <>
              <div style={{
                display: 'grid',
                gridTemplateColumns: 'repeat(2, 1fr)',
                gap: '1rem',
                marginBottom: '1rem'
              }}>
                <div style={{
                  background: 'var(--bg-secondary)',
                  padding: '1rem',
                  borderRadius: '8px',
                  textAlign: 'center'
                }}>
                  <div style={{
                    fontSize: '1.75rem',
                    fontWeight: 700,
                    background: 'linear-gradient(135deg, var(--accent-primary), var(--accent-secondary))',
                    WebkitBackgroundClip: 'text',
                    WebkitTextFillColor: 'transparent'
                  }}>
                    {result.qc_results?.data?.n_cells?.toLocaleString() || 'N/A'}
                  </div>
                  <div style={{ color: 'var(--text-secondary)', fontSize: '0.9rem' }}>Cells</div>
                </div>
                <div style={{
                  background: 'var(--bg-secondary)',
                  padding: '1rem',
                  borderRadius: '8px',
                  textAlign: 'center'
                }}>
                  <div style={{
                    fontSize: '1.75rem',
                    fontWeight: 700,
                    background: 'linear-gradient(135deg, var(--accent-primary), var(--accent-secondary))',
                    WebkitBackgroundClip: 'text',
                    WebkitTextFillColor: 'transparent'
                  }}>
                    {result.qc_results?.data?.n_genes?.toLocaleString() || 'N/A'}
                  </div>
                  <div style={{ color: 'var(--text-secondary)', fontSize: '0.9rem' }}>Genes</div>
                </div>
              </div>

              {result.validation && (
                <div style={{
                  padding: '0.75rem 1rem',
                  borderRadius: '8px',
                  background: result.validation.decision === 'PROCEED'
                    ? 'rgba(16, 185, 129, 0.1)'
                    : 'rgba(251, 191, 36, 0.1)',
                  border: `1px solid ${result.validation.decision === 'PROCEED' ? 'var(--success)' : '#fbbf24'}`,
                  marginBottom: '1rem'
                }}>
                  <strong>Validation: </strong>
                  {result.validation.decision === 'PROCEED' && '✅ All checks passed'}
                  {result.validation.decision === 'PROCEED_WITH_CAUTION' && '⚠️ Proceed with caution'}
                  {result.validation.warnings?.length > 0 && (
                    <ul style={{ margin: '0.5rem 0 0 1rem', fontSize: '0.9rem' }}>
                      {result.validation.warnings.map((w, i) => <li key={i}>{w}</li>)}
                    </ul>
                  )}
                </div>
              )}

              <p style={{ fontSize: '0.9rem', color: 'var(--text-secondary)' }}>
                Navigate to the <strong>Pipeline</strong> tab to configure and run the full analysis.
              </p>
            </>
          )}

          <style>{`
            @keyframes fadeIn {
              from { opacity: 0; transform: translateY(10px); }
              to { opacity: 1; transform: translateY(0); }
            }
          `}</style>
        </div>
      )}

      {status === 'error' && (
        <div style={{
          color: 'var(--error)',
          background: 'rgba(239, 68, 68, 0.1)',
          padding: '1rem',
          borderRadius: '8px'
        }}>
          <strong>Error:</strong> {error}
        </div>
      )}
    </div>
  );
}
