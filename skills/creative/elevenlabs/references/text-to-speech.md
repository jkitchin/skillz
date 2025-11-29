

# ElevenLabs Text-to-Speech Reference

## Overview

ElevenLabs Text-to-Speech (TTS) converts written text into lifelike audio with natural intonation, emotional awareness, and support for 32 languages. The API provides three main models optimized for different use cases, 100+ professional voices, and advanced features like streaming and multi-speaker generation.

## Available Models

### Eleven Multilingual v2

**Model ID:** `eleven_multilingual_v2`

**Specifications:**
- Languages: 29 languages
- Quality: Highest available
- Latency: Standard (higher than Flash)
- Use Case: Final deliverables, audiobooks, professional content

**Strengths:**
- Best audio quality and naturalness
- Most nuanced emotional expression
- Excellent for long-form content
- Superior pronunciation accuracy

**Best For:**
- Audiobooks and narration
- Professional voice-overs
- Content where quality is paramount
- Characters and storytelling

### Eleven Flash v2.5

**Model ID:** `eleven_flash_v2_5`

**Specifications:**
- Languages: 32 languages
- Quality: High (slightly below Multilingual v2)
- Latency: Ultra-low (75ms)
- Cost: 50% lower than Multilingual v2
- Use Case: Real-time applications, streaming

**Strengths:**
- Extremely low latency
- Cost-effective
- Suitable for real-time generation
- Good quality/speed balance

**Best For:**
- Conversational AI
- Live streaming applications
- Real-time voice synthesis
- High-volume use cases

### Eleven Turbo v2.5

**Model ID:** `eleven_turbo_v2_5`

**Specifications:**
- Languages: 32 languages
- Quality: Balanced
- Latency: Low (between Flash and Multilingual)
- Use Case: General developer applications

**Strengths:**
- Balanced latency and quality
- Good for most applications
- Versatile across use cases

**Best For:**
- General applications
- Video content
- Podcasts
- E-learning materials

## API Usage

### Basic Setup

```python
import os
from elevenlabs.client import ElevenLabs
from pathlib import Path

# Initialize client
client = ElevenLabs(api_key=os.environ.get("ELEVENLABS_API_KEY"))
```

### Simple Text-to-Speech

```python
# Convert text to speech
audio = client.text_to_speech.convert(
    text="Hello, this is a test of the text to speech system.",
    voice_id="JBFqnCBsd6RMkjVDRZzb",  # George voice
    model_id="eleven_multilingual_v2",
    output_format="mp3_44100_128"
)

# Save to file
output_path = Path("output.mp3")
with output_path.open("wb") as f:
    for chunk in audio:
        f.write(chunk)
```

### Streaming Audio

```python
from elevenlabs import stream

# Generate and stream audio in real-time
audio_stream = client.text_to_speech.convert_as_stream(
    text="This audio is being streamed as it generates.",
    voice_id="JBFqnCBsd6RMkjVDRZzb",
    model_id="eleven_flash_v2_5",  # Use Flash for low latency
    output_format="mp3_44100_128"
)

# Stream to speakers
stream(audio_stream)
```

### Save Streamed Audio

```python
# Stream and save simultaneously
audio_stream = client.text_to_speech.convert_as_stream(
    text="Long text content here...",
    voice_id="21m00Tcm4TlvDq8ikWAM",  # Rachel
    model_id="eleven_flash_v2_5"
)

output_path = Path("streamed_output.mp3")
with output_path.open("wb") as f:
    for chunk in audio_stream:
        f.write(chunk)
```

## Voice Management

### List All Voices

```python
# Get all available voices
voices_response = client.voices.get_all()

for voice in voices_response.voices:
    print(f"Name: {voice.name}")
    print(f"ID: {voice.voice_id}")
    print(f"Category: {voice.category}")
    print(f"Description: {voice.description}")
    print(f"Labels: {voice.labels}")
    print("---")
```

### Search Voices by Criteria

```python
# Filter voices by characteristics
voices_response = client.voices.get_all()

# Find female voices
female_voices = [v for v in voices_response.voices
                  if 'female' in str(v.labels).lower()]

# Find young voices
young_voices = [v for v in voices_response.voices
                 if 'young' in str(v.labels).lower()]

# Find voices by language
english_voices = [v for v in voices_response.voices
                   if 'english' in str(v.labels).lower()]
```

### Get Voice Settings

```python
# Get settings for a specific voice
voice_settings = client.voices.get_settings(
    voice_id="JBFqnCBsd6RMkjVDRZzb"
)

print(f"Stability: {voice_settings.stability}")
print(f"Similarity Boost: {voice_settings.similarity_boost}")
```

## Common Voice IDs

### English Voices

**Male Voices:**
- `JBFqnCBsd6RMkjVDRZzb` - **George** (middle-aged, narrative)
- `TxGEqnHWrfWFTfGW9XjX` - **Josh** (young, energetic)
- `ErXwobaYiN019PkySvjV` - **Antoni** (young, well-rounded)
- `VR6AewLTigWG4xSOukaG` - **Arnold** (deep, mature)
- `pNInz6obpgDQGcFmaJgB` - **Adam** (deep, authoritative)

**Female Voices:**
- `21m00Tcm4TlvDq8ikWAM` - **Rachel** (young, calm)
- `AZnzlk1XvdvUeBnXmlld` - **Domi** (young, confident)
- `EXAVITQu4vr4xnSDxMaL` - **Bella** (young, expressive)
- `MF3mGyEYCl7XYWbV9V6O` - **Elli** (young, emotional)
- `XrExE9yKIg1WjnnlVkGX` - **Matilda** (young, warm)

### Multilingual Voices

- `ThT5KcBeYPX3keUQqHPh` - **Freya** (Spanish, female)
- `XB0fDUnXU5powFXDhCwa` - **Charlotte** (Swedish, female)
- `iP95p4xoKVk53GoZ742B` - **Chris** (American, male)

## Output Formats

### Format Specification

Format string pattern: `{codec}_{sample_rate}_{bitrate}`

### MP3 Formats

**High Quality:**
- `mp3_44100_192` - 44.1kHz, 192kbps (highest quality)
- `mp3_44100_128` - 44.1kHz, 128kbps (standard high quality)

**Standard Quality:**
- `mp3_44100_96` - 44.1kHz, 96kbps
- `mp3_44100_64` - 44.1kHz, 64kbps

**Lower Quality (smaller files):**
- `mp3_22050_32` - 22.05kHz, 32kbps

### PCM Formats

**CD Quality:**
- `pcm_44100` - 44.1kHz, 16-bit (lossless, large files)

**Other Sample Rates:**
- `pcm_24000` - 24kHz
- `pcm_22050` - 22.05kHz
- `pcm_16000` - 16kHz
- `pcm_8000` - 8kHz (telephony)

### Other Formats

**Opus:**
- `opus_48000_192` - 48kHz, 192kbps
- `opus_48000_128` - 48kHz, 128kbps

**Telephony:**
- `ulaw_8000` - 8kHz µ-law
- `alaw_8000` - 8kHz A-law

### Format Selection Guide

- **Web/Mobile:** `mp3_44100_128` (good balance)
- **High Quality:** `mp3_44100_192` or `pcm_44100`
- **Streaming:** `mp3_44100_96` (smaller chunks)
- **Voice Calls:** `ulaw_8000` or `pcm_8000`
- **Maximum Quality:** `pcm_44100` (lossless)

## Advanced Features

### Multi-Speaker Dialogue

```python
from pydub import AudioSegment

# Define conversation
dialogue = [
    ("JBFqnCBsd6RMkjVDRZzb", "Hello, how can I help you today?"),
    ("21m00Tcm4TlvDq8ikWAM", "I'm looking for information about your products."),
    ("JBFqnCBsd6RMkjVDRZzb", "I'd be happy to help with that!"),
]

# Generate each line
combined = AudioSegment.empty()

for voice_id, text in dialogue:
    audio = client.text_to_speech.convert(
        text=text,
        voice_id=voice_id,
        model_id="eleven_multilingual_v2"
    )

    # Save temporarily
    temp_file = Path(f"temp_{voice_id}.mp3")
    with temp_file.open("wb") as f:
        for chunk in audio:
            f.write(chunk)

    # Add to combined audio
    segment = AudioSegment.from_mp3(str(temp_file))
    combined += segment

    # Add pause between speakers (500ms)
    combined += AudioSegment.silent(duration=500)

    # Clean up
    temp_file.unlink()

# Export final dialogue
combined.export("conversation.mp3", format="mp3")
```

### Long-Form Content

```python
def split_text(text, max_length=5000):
    """Split long text into chunks at sentence boundaries"""
    sentences = text.split('. ')
    chunks = []
    current_chunk = ""

    for sentence in sentences:
        if len(current_chunk) + len(sentence) < max_length:
            current_chunk += sentence + ". "
        else:
            chunks.append(current_chunk.strip())
            current_chunk = sentence + ". "

    if current_chunk:
        chunks.append(current_chunk.strip())

    return chunks

# Process long text
long_text = "... very long text here ..."
chunks = split_text(long_text)

combined = AudioSegment.empty()

for i, chunk in enumerate(chunks):
    print(f"Processing chunk {i+1}/{len(chunks)}...")

    audio = client.text_to_speech.convert(
        text=chunk,
        voice_id="JBFqnCBsd6RMkjVDRZzb",
        model_id="eleven_multilingual_v2"
    )

    # Save chunk temporarily
    temp_file = Path(f"chunk_{i}.mp3")
    with temp_file.open("wb") as f:
        for data in audio:
            f.write(data)

    # Add to combined
    segment = AudioSegment.from_mp3(str(temp_file))
    combined += segment
    temp_file.unlink()

# Export complete audio
combined.export("audiobook.mp3", format="mp3")
```

### Async Generation

```python
from elevenlabs.client import AsyncElevenLabs
import asyncio

async def generate_multiple_audio():
    """Generate multiple audio files concurrently"""
    client = AsyncElevenLabs(api_key=os.environ["ELEVENLABS_API_KEY"])

    texts = [
        "First audio clip",
        "Second audio clip",
        "Third audio clip"
    ]

    tasks = []
    for i, text in enumerate(texts):
        task = client.text_to_speech.convert(
            text=text,
            voice_id="JBFqnCBsd6RMkjVDRZzb",
            model_id="eleven_flash_v2_5"
        )
        tasks.append((i, task))

    # Execute concurrently
    for i, task in tasks:
        audio = await task
        output_path = Path(f"audio_{i}.mp3")
        with output_path.open("wb") as f:
            async for chunk in audio:
                f.write(chunk)

# Run async function
asyncio.run(generate_multiple_audio())
```

## Voice Settings & Customization

### Voice Parameters

**Stability (0.0 - 1.0):**
- Low (0.0-0.3): More variable, expressive, emotional
- Medium (0.4-0.6): Balanced consistency and expression
- High (0.7-1.0): More consistent, stable, monotone

**Similarity Boost (0.0 - 1.0):**
- Low (0.0-0.3): More creative interpretation
- Medium (0.4-0.6): Balanced
- High (0.7-1.0): Closer to original voice sample

**Style (0.0 - 1.0):**
- Controls style exaggeration
- Higher values = more dramatic delivery

**Use Speaker Boost (boolean):**
- Enhances voice similarity
- Useful for voice cloning

### Custom Settings Example

```python
# Generate with custom voice settings
audio = client.text_to_speech.convert(
    text="This uses custom voice settings.",
    voice_id="JBFqnCBsd6RMkjVDRZzb",
    model_id="eleven_multilingual_v2",
    voice_settings={
        "stability": 0.5,
        "similarity_boost": 0.75,
        "style": 0.0,
        "use_speaker_boost": True
    }
)
```

## Supported Languages

ElevenLabs supports 32 languages across models:

**Western European:**
English, Spanish, French, German, Italian, Portuguese, Dutch, Polish, Swedish

**Eastern European:**
Czech, Russian, Ukrainian, Romanian, Bulgarian, Croatian, Slovak

**Asian:**
Chinese (Mandarin), Japanese, Korean, Hindi, Tamil, Telugu, Filipino, Malay, Indonesian, Vietnamese, Thai

**Other:**
Arabic, Turkish

**Usage:**
```python
# Generate in Spanish
audio = client.text_to_speech.convert(
    text="Hola, ¿cómo estás?",
    voice_id="ThT5KcBeYPX3keUQqHPh",  # Spanish voice
    model_id="eleven_multilingual_v2"
)
```

## Best Practices

### Model Selection

1. **Audiobooks/Narration** → Multilingual v2 (highest quality)
2. **Real-time Apps** → Flash v2.5 (low latency)
3. **General Use** → Turbo v2.5 (balanced)
4. **Cost-Sensitive** → Flash v2.5 (50% cheaper)

### Voice Selection

1. **Match content**: Choose voice appropriate for audience and tone
2. **Test multiple**: Generate samples with different voices
3. **Consider labels**: Use voice labels (age, accent, style) to filter
4. **Consistency**: Use same voice for series/brand content

### Quality Optimization

1. **Clean text**: Remove formatting artifacts, fix typos
2. **Proper punctuation**: Affects phrasing and pauses
3. **SSML tags**: Use for fine control (future feature)
4. **Split long texts**: Avoid generation timeouts
5. **High bitrate**: Use 128kbps+ for professional content

### Performance

1. **Stream for UI**: Use streaming for better user experience
2. **Batch processing**: Generate multiple short clips concurrently
3. **Cache audio**: Reuse generated audio when possible
4. **Monitor quota**: Track character usage to avoid limits

### Error Handling

```python
import time

def generate_with_retry(client, text, voice_id, max_retries=3):
    """Generate audio with retry logic"""
    for attempt in range(max_retries):
        try:
            return client.text_to_speech.convert(
                text=text,
                voice_id=voice_id,
                model_id="eleven_multilingual_v2"
            )
        except Exception as e:
            if "rate limit" in str(e).lower() and attempt < max_retries - 1:
                wait_time = 2 ** attempt  # Exponential backoff
                print(f"Rate limited. Waiting {wait_time}s...")
                time.sleep(wait_time)
            else:
                raise
```

## Troubleshooting

### Common Issues

**Poor Audio Quality:**
- Use higher bitrate format (`mp3_44100_128` or higher)
- Try different voice
- Use Multilingual v2 model
- Check text for formatting issues

**Unnatural Speech:**
- Adjust voice settings (lower stability for more expression)
- Fix text punctuation and capitalization
- Split very long sentences
- Use appropriate voice for content type

**Generation Errors:**
- Verify API key is set correctly
- Check quota/usage limits
- Ensure text is within character limits
- Remove special characters if causing issues

**Slow Generation:**
- Use Flash v2.5 for faster results
- Consider streaming for immediate feedback
- Check internet connection
- Process in chunks for long content

## API Limits

### Free Tier
- 10,000 characters per month
- Access to all voices
- Text-to-speech only (no music)

### Paid Tiers
- Higher character limits
- Priority processing
- Music generation access
- Commercial usage rights

Check current limits at: https://elevenlabs.io/pricing

## Resources

- **API Documentation**: https://elevenlabs.io/docs/api-reference/text-to-speech
- **Voice Library**: https://elevenlabs.io/voice-library
- **Usage Dashboard**: https://elevenlabs.io/app/usage
- **Python SDK**: https://github.com/elevenlabs/elevenlabs-python
- **Supported Languages**: https://elevenlabs.io/languages
