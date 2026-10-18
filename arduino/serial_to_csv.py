"""Save the Arduino logger's serial output to a CSV file.

Usage:
    pip install pyserial
    python serial_to_csv.py --port COM3 --out data/raw/glycerol_10wt_run1.csv
    (on Linux/macOS the port looks like /dev/ttyUSB0 or /dev/cu.usbmodem...)

Stop with Ctrl+C. Each line is also echoed to the screen so you can watch the run.
Name files as  <solute>_<wt>wt_run<n>.csv  so the analysis can parse them.
"""

import argparse
import sys
import time

try:
    import serial
except ImportError:
    sys.exit("pyserial not installed: pip install pyserial")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--port", required=True)
    ap.add_argument("--baud", type=int, default=115200)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    with serial.Serial(a.port, a.baud, timeout=5) as ser, open(a.out, "w") as f:
        time.sleep(2)                      # Arduino resets when the port opens
        ser.reset_input_buffer()
        print(f"Logging {a.port} -> {a.out}  (Ctrl+C to stop)")
        try:
            while True:
                line = ser.readline().decode(errors="ignore").strip()
                if not line:
                    continue
                print(line)
                f.write(line + "\n")
                f.flush()
        except KeyboardInterrupt:
            print("\nStopped.")


if __name__ == "__main__":
    main()
