function initializeOscillatingField() {
    const fieldCanvas = document.getElementById('fieldCanvas');
    const graphCtx = fieldCanvas.getContext('2d');

    // Hardcoded parameters
    const G = 1;
    const H = 0;
    const A = 1;
    const omega = 2;
    const phi = 0;
    const sigma = 15;
    const t0 = 65;
    const normalizationTerm = 1; // Placeholder, adjust as needed
    const maxTime = 125;
    const timeStep = 0.05;

    const width = fieldCanvas.width;
    const height = fieldCanvas.height;

    let time = 0;
    const EValues = [];
    const timeValues = [];

    function gaussianPulse(t, t0, sigma, G, H) {
        const tDiff = t - t0;
        const gaussian = (G / normalizationTerm) * Math.exp(-Math.pow(tDiff, 2) / (2 * Math.pow(sigma, 2)));
        const heaviside = H * (tDiff > 0 ? 1 : 0);
        return gaussian + heaviside;
    }

    function oscillation(t, t0, A, omega, phi, C = 0) {
        return A * Math.sin(omega * (t - t0) - phi) + C;
    }

    function calculateElectricField(t) {
        const term1 = gaussianPulse(t, t0, sigma, G, H);
        const term2 = oscillation(t, t0, A, omega, phi);
        return term1 * term2;
    }


    function drawField(currentTime, currentE) {
        graphCtx.clearRect(0, 0, width, height);

        graphCtx.beginPath();
        graphCtx.moveTo(0, height / 2);

        for (let i = 0; i < EValues.length; i++) {
            const x = (i / EValues.length) * width;
            const y = height / 2 - EValues[i] * 100; // Scale for visibility
            graphCtx.lineTo(x, y);
        }

        graphCtx.strokeStyle = 'red';
        graphCtx.lineWidth = 2;
        graphCtx.stroke();

        // Add current time and field value
        graphCtx.font = '16px Arial';
        graphCtx.fillStyle = 'black';
        graphCtx.fillText(`Time: ${currentTime.toFixed(2)} fs`, 10, 20);
        graphCtx.fillText(`Field Value: ${currentE.toFixed(3)}`, 10, 40);
    }

    function animate() {
        if (time >= maxTime) return;

        const E = calculateElectricField(time);
        EValues.push(E);
        timeValues.push(time);

        drawField(time, E);

        time += timeStep;
        requestAnimationFrame(animate);
    }

    function restartAnimation() {
        time = 0;
        EValues.length = 0;
        timeValues.length = 0;
        animate();
    }

    document.getElementById('restartButton').addEventListener('click', restartAnimation);

    restartAnimation();
}
initializeOscillatingField();
