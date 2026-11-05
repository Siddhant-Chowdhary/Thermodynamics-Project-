// Cooling-curve logger: up to 3 DS18B20 probes on one data pin.
// Prints CSV to serial:  time_s,T1,T2,T3   every LOG_INTERVAL_MS.
//
// Wiring (each probe): red -> 5V, black -> GND, yellow (data) -> pin 2,
// with ONE 4.7 kOhm pull-up resistor between data and 5V (shared by all probes).
//
// Libraries (Arduino IDE > Library Manager): "OneWire" and "DallasTemperature".

#include <OneWire.h>
#include <DallasTemperature.h>

const int ONE_WIRE_PIN = 2;
const unsigned long LOG_INTERVAL_MS = 2000;   // 2 s between readings
const int MAX_PROBES = 3;

OneWire oneWire(ONE_WIRE_PIN);
DallasTemperature sensors(&oneWire);
DeviceAddress addr[MAX_PROBES];
int nProbes = 0;
unsigned long t0;

void setup() {
  Serial.begin(115200);
  sensors.begin();
  nProbes = min((int)sensors.getDeviceCount(), MAX_PROBES);
  for (int i = 0; i < nProbes; i++) {
    sensors.getAddress(addr[i], i);
    sensors.setResolution(addr[i], 12);        // 0.0625 C steps, ~750 ms conversion
  }
  sensors.setWaitForConversion(true);
  Serial.print("time_s");
  for (int i = 0; i < nProbes; i++) { Serial.print(",T"); Serial.print(i + 1); }
  Serial.println();
  t0 = millis();
}

void loop() {
  static unsigned long next = 0;
  unsigned long now = millis();
  if (now < next) return;
  next = now + LOG_INTERVAL_MS;

  sensors.requestTemperatures();
  Serial.print((now - t0) / 1000.0, 1);
  for (int i = 0; i < nProbes; i++) {
    Serial.print(",");
    float c = sensors.getTempC(addr[i]);
    if (c == DEVICE_DISCONNECTED_C) Serial.print("nan");
    else Serial.print(c, 4);
  }
  Serial.println();
}
