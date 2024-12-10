const int inputPin = 12;
const int outputPin = 13;
const int time = 100;

void setup() {
  pinMode(inputPin, INPUT);
  pinMode(outputPin, OUTPUT);

}

void loop() {
  bool inputState = digitalRead(inputPin) == HIGH;
  if (inputState) {
    digitalWrite(outputPin, HIGH);
    delay(time);
    digitalWrite(outputPin, LOW);
    delay(1000-time);
  }
}
