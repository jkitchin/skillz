# ElevenLabs Complete Examples

## Overview

This guide provides complete, runnable examples for common ElevenLabs use cases covering text-to-speech, sound effects, and music generation.

## Setup

```python
import os
from elevenlabs.client import ElevenLabs
from pathlib import Path

# Initialize client
client = ElevenLabs(api_key=os.environ.get("ELEVENLABS_API_KEY"))

# Helper function to save audio
def save_audio(audio_generator, filename):
    """Save audio generator to file"""
    output_path = Path(filename)
    with output_path.open("wb") as f:
        for chunk in audio_generator:
            f.write(chunk)
    print(f"Saved: {output_path.absolute()}")
    return output_path
```

## Example 1: Audiobook Narration

**Use Case:** Convert a chapter of text into professional audiobook format

```python
# Chapter text
chapter_text = """
Chapter One: The Beginning

It was a dark and stormy night. The rain fell in sheets against the windows
of the old mansion, and thunder rolled across the distant hills. Inside,
a single candle flickered in the library, casting long shadows across the
rows of ancient books.

Sarah clutched her grandmother's letter tightly as she made her way down
the creaking hallway. The words echoed in her mind: "When you read this,
I will be gone. But the secrets of Blackwood Manor must not die with me."
"""

# Generate audiobook narration
audio = client.text_to_speech.convert(
    text=chapter_text,
    voice_id="JBFqnCBsd6RMkjVDRZzb",  # George - narrative voice
    model_id="eleven_multilingual_v2",  # Highest quality
    output_format="mp3_44100_128"
)

# Save chapter
save_audio(audio, "chapter_01.mp3")
print("Audiobook chapter generated successfully!")
```

## Example 2: Multi-Language Product Demo

**Use Case:** Create product welcome messages in multiple languages

```python
# Welcome messages in different languages
messages = {
    "english": {
        "text": "Welcome to our application! We're excited to have you here.",
        "voice_id": "21m00Tcm4TlvDq8ikWAM",  # Rachel
    },
    "spanish": {
        "text": "¡Bienvenido a nuestra aplicación! Estamos emocionados de tenerte aquí.",
        "voice_id": "ThT5KcBeYPX3keUQqHPh",  # Spanish voice
    },
    "french": {
        "text": "Bienvenue dans notre application! Nous sommes ravis de vous accueillir.",
        "voice_id": "XB0fDUnXU5powFXDhCwa",  # French voice
    }
}

# Generate each language version
for lang, config in messages.items():
    audio = client.text_to_speech.convert(
        text=config["text"],
        voice_id=config["voice_id"],
        model_id="eleven_multilingual_v2"
    )
    save_audio(audio, f"welcome_{lang}.mp3")

print("Multi-language welcome messages generated!")
```

## Example 3: Video Game Sound Effects Library

**Use Case:** Generate complete sound effects library for a fantasy RPG

```python
# Define game sound effects
game_sfx = {
    "combat": [
        ("sword_swing", "sword whooshing through air, fast swing, combat sound", 1.0),
        ("shield_block", "metal shield blocking sword attack, metallic clang", 0.8),
        ("bow_shot", "arrow being shot from bow, twang and whoosh", 1.2),
        ("spell_fireball", "magical fireball spell launching, fire whoosh", 1.5),
    ],
    "ui": [
        ("menu_open", "fantasy game menu opening, magical chime", 0.6),
        ("item_pickup", "picking up item, positive collect sound", 0.4),
        ("quest_complete", "quest completion fanfare, triumphant", 2.0),
        ("level_up", "character level up sound, ascending magical tones", 2.5),
    ],
    "environment": [
        ("footsteps_stone", "footsteps on stone dungeon floor, echoing", 3.0),
        ("door_creak", "old wooden door creaking open, dungeon atmosphere", 2.5),
        ("torch_fire", "burning torch crackling, ambient fire sound", 5.0),
        ("wind_howl", "wind howling through castle, eerie atmosphere", 8.0),
    ]
}

# Generate all sound effects organized by category
for category, effects in game_sfx.items():
    category_dir = Path(f"game_sfx/{category}")
    category_dir.mkdir(parents=True, exist_ok=True)

    print(f"\nGenerating {category} sounds...")
    for name, description, duration in effects:
        audio = client.text_to_sound_effects.convert(
            text=description,
            duration_seconds=duration
        )
        save_audio(audio, category_dir / f"{name}.mp3")

print("\nGame sound effects library generated!")
```

## Example 4: Podcast Production Suite

**Use Case:** Create complete audio package for a podcast episode

```python
from pydub import AudioSegment

# 1. Generate intro music
print("Generating intro music...")
intro_music = client.music_generation.compose(
    prompt="""Upbeat podcast intro music with modern electronic beats,
    catchy synth melody, and professional sound. Energetic and memorable.""",
    music_length_ms=15000  # 15 seconds
)
save_audio(intro_music, "podcast_intro_music.mp3")

# 2. Generate intro voiceover
print("Generating intro voiceover...")
intro_vo = client.text_to_speech.convert(
    text="""Welcome to Tech Talk, the podcast where we dive deep into the
    latest technology trends and innovations!""",
    voice_id="TxGEqnHWrfWFTfGW9XjX",  # Josh - energetic
    model_id="eleven_flash_v2_5"
)
save_audio(intro_vo, "podcast_intro_voice.mp3")

# 3. Generate transition sounds
print("Generating transitions...")
transitions = [
    ("segment_transition", "smooth podcast transition sound, modern swoosh", 2.0),
    ("ad_break_in", "ad break intro, upbeat jingle", 3.0),
    ("ad_break_out", "ad break outro, return to show", 3.0),
]

for name, description, duration in transitions:
    audio = client.text_to_sound_effects.convert(
        text=description,
        duration_seconds=duration
    )
    save_audio(audio, f"podcast_{name}.mp3")

# 4. Generate outro
print("Generating outro...")
outro_vo = client.text_to_speech.convert(
    text="""Thanks for listening to Tech Talk! Don't forget to subscribe
    and join us next week for more tech insights!""",
    voice_id="TxGEqnHWrfWFTfGW9XjX",
    model_id="eleven_flash_v2_5"
)
save_audio(outro_vo, "podcast_outro_voice.mp3")

print("\nPodcast production suite complete!")
print("Note: Use audio editing software to mix voice and music")
```

## Example 5: E-Learning Course Narration

**Use Case:** Generate narration for online course modules with multiple speakers

```python
# Course content with dialogue
course_content = [
    {
        "speaker": "instructor",
        "voice_id": "JBFqnCBsd6RMkjVDRZzb",  # George
        "text": """Welcome to Module 3: Python Basics. In this module,
        we'll learn about variables, data types, and basic operations."""
    },
    {
        "speaker": "student",
        "voice_id": "21m00Tcm4TlvDq8ikWAM",  # Rachel
        "text": """That sounds interesting! Can you give us an example?"""
    },
    {
        "speaker": "instructor",
        "voice_id": "JBFqnCBsd6RMkjVDRZzb",
        "text": """Of course! Let's start with a simple variable assignment.
        In Python, you can create a variable by simply assigning a value to a name."""
    },
]

# Generate each segment
from pydub import AudioSegment

combined = AudioSegment.empty()

for i, segment in enumerate(course_content):
    print(f"Generating segment {i+1}/{len(course_content)}...")

    audio = client.text_to_speech.convert(
        text=segment["text"],
        voice_id=segment["voice_id"],
        model_id="eleven_multilingual_v2"
    )

    # Save temporarily
    temp_file = Path(f"temp_{i}.mp3")
    save_audio(audio, temp_file)

    # Add to combined audio
    audio_segment = AudioSegment.from_mp3(str(temp_file))
    combined += audio_segment

    # Add pause between speakers
    combined += AudioSegment.silent(duration=800)

    # Clean up temp file
    temp_file.unlink()

# Export final course audio
combined.export("course_module_3.mp3", format="mp3")
print("Course module narration complete!")
```

## Example 6: YouTube Video Background Package

**Use Case:** Create complete audio package for a YouTube video

```python
# 1. Generate background music
print("Generating background music...")
bg_music = client.music_generation.compose(
    prompt="""Calm tutorial background music with soft piano and ambient pads,
    non-distracting and focused, perfect for educational content.""",
    music_length_ms=180000  # 3 minutes
)
save_audio(bg_music, "video_background_music.mp3")

# 2. Generate intro/outro sounds
print("Generating intro/outro...")
sounds = [
    ("intro_whoosh", "YouTube video intro whoosh, energetic and modern", 2.0),
    ("subscribe_bell", "Subscribe notification bell sound, clear and pleasant", 1.0),
    ("outro_whoosh", "Video outro transition, smooth ending", 2.5),
]

for name, description, duration in sounds:
    audio = client.text_to_sound_effects.convert(
        text=description,
        duration_seconds=duration
    )
    save_audio(audio, f"video_{name}.mp3")

# 3. Generate voiceover
print("Generating voiceover...")
script = """Hey everyone, welcome back to the channel! Today we're going to
learn about Python programming. If you're new here, don't forget to hit that
subscribe button and turn on notifications!"""

voiceover = client.text_to_speech.convert(
    text=script,
    voice_id="TxGEqnHWrfWFTfGW9XjX",  # Josh - enthusiastic
    model_id="eleven_multilingual_v2"
)
save_audio(voiceover, "video_voiceover.mp3")

print("\nYouTube video audio package complete!")
```

## Example 7: Mobile App Sound Design

**Use Case:** Generate UI sounds for a mobile application

```python
# Define app UI sounds
ui_sounds = {
    "navigation": [
        ("button_tap", "soft UI button tap, satisfying click", 0.2),
        ("swipe", "smooth swipe gesture sound, fluid transition", 0.4),
        ("back_button", "back navigation sound, subtle pop", 0.3),
    ],
    "feedback": [
        ("success", "success action sound, positive chime", 0.6),
        ("error", "error notification, gentle warning tone", 0.7),
        ("notification", "new notification sound, attention-grabbing but pleasant", 1.0),
    ],
    "actions": [
        ("send_message", "message sent sound, whoosh and confirm", 0.8),
        ("receive_message", "message received, subtle ping", 0.5),
        ("pull_refresh", "pull to refresh sound, water drop effect", 1.0),
        ("upload_complete", "file upload complete, ascending tones", 1.5),
    ]
}

# Generate all UI sounds
for category, sounds in ui_sounds.items():
    category_dir = Path(f"app_sounds/{category}")
    category_dir.mkdir(parents=True, exist_ok=True)

    print(f"\nGenerating {category} sounds...")
    for name, description, duration in sounds:
        audio = client.text_to_sound_effects.convert(
            text=description,
            duration_seconds=duration,
            prompt_influence=0.6  # Balanced control
        )
        save_audio(audio, category_dir / f"{name}.mp3")

print("\nMobile app sound design complete!")
```

## Example 8: Meditation/Relaxation App

**Use Case:** Create audio content for wellness application

```python
# 1. Generate ambient background music
print("Generating ambient music...")
ambient_tracks = [
    ("forest", """Peaceful forest ambience music with gentle bird sounds,
    flowing water, and soft wind chimes. Calming and natural.""", 300000),

    ("ocean", """Ocean waves meditation music with gentle waves, soft pads,
    and peaceful atmosphere. Relaxing and soothing.""", 300000),

    ("rain", """Rain meditation music with gentle rain sounds, distant thunder,
    and calming piano. Peaceful sleep music.""", 300000),
]

for name, prompt, duration in ambient_tracks:
    audio = client.music_generation.compose(
        prompt=prompt,
        music_length_ms=duration
    )
    save_audio(audio, f"meditation_{name}.mp3")

# 2. Generate guided meditation narration
print("\nGenerating guided meditation...")
meditation_script = """
Find a comfortable position and gently close your eyes. Take a deep breath in
through your nose... and slowly exhale through your mouth. With each breath,
feel your body becoming more relaxed. Let go of any tension in your shoulders...
your neck... your jaw. You are safe. You are calm. You are at peace.
"""

meditation_voice = client.text_to_speech.convert(
    text=meditation_script,
    voice_id="21m00Tcm4TlvDq8ikWAM",  # Rachel - calm voice
    model_id="eleven_multilingual_v2",
    voice_settings={
        "stability": 0.7,  # Very stable, calm delivery
        "similarity_boost": 0.75
    }
)
save_audio(meditation_voice, "guided_meditation.mp3")

# 3. Generate sound effects
print("\nGenerating nature sounds...")
nature_sounds = [
    ("singing_bowl", "Tibetan singing bowl sound, meditation bell, resonant", 5.0),
    ("chimes", "Wind chimes tinkling gently, peaceful garden sounds", 8.0),
    ("brook", "Gentle brook babbling over rocks, forest stream", 10.0),
]

for name, description, duration in nature_sounds:
    audio = client.text_to_sound_effects.convert(
        text=description,
        duration_seconds=duration
    )
    save_audio(audio, f"nature_{name}.mp3")

print("\nMeditation app audio content complete!")
```

## Example 9: Batch Processing for Content Library

**Use Case:** Generate large library of content efficiently

```python
import time

def batch_generate_with_retry(items, item_type="speech"):
    """
    Generate multiple items with error handling and retry logic
    """
    results = {"success": [], "failed": []}

    for i, item in enumerate(items, 1):
        print(f"\n[{i}/{len(items)}] Processing: {item['name']}")

        max_retries = 3
        for attempt in range(max_retries):
            try:
                if item_type == "speech":
                    audio = client.text_to_speech.convert(
                        text=item["text"],
                        voice_id=item.get("voice_id", "JBFqnCBsd6RMkjVDRZzb"),
                        model_id=item.get("model", "eleven_multilingual_v2")
                    )
                elif item_type == "sound":
                    audio = client.text_to_sound_effects.convert(
                        text=item["description"],
                        duration_seconds=item.get("duration", 3.0)
                    )
                elif item_type == "music":
                    audio = client.music_generation.compose(
                        prompt=item["prompt"],
                        music_length_ms=item.get("duration_ms", 30000)
                    )

                # Save
                output_path = save_audio(audio, item["filename"])
                results["success"].append((item["name"], output_path))
                break  # Success, exit retry loop

            except Exception as e:
                if "rate limit" in str(e).lower() and attempt < max_retries - 1:
                    wait_time = 2 ** attempt
                    print(f"  Rate limited. Waiting {wait_time}s...")
                    time.sleep(wait_time)
                elif attempt == max_retries - 1:
                    results["failed"].append((item["name"], str(e)))
                    print(f"  Failed after {max_retries} attempts: {e}")
                else:
                    results["failed"].append((item["name"], str(e)))
                    break

        # Small delay between items
        time.sleep(0.5)

    return results

# Example: Batch generate speech for app notifications
notifications = [
    {
        "name": "welcome",
        "text": "Welcome! Your account has been created successfully.",
        "filename": "notif_welcome.mp3"
    },
    {
        "name": "order_shipped",
        "text": "Great news! Your order has been shipped.",
        "filename": "notif_shipped.mp3"
    },
    {
        "name": "reminder",
        "text": "Don't forget! You have a meeting in 15 minutes.",
        "filename": "notif_reminder.mp3"
    },
]

results = batch_generate_with_retry(notifications, item_type="speech")

print(f"\n✓ Successfully generated: {len(results['success'])}")
print(f"✗ Failed: {len(results['failed'])}")
```

## Example 10: Complete Video Production Workflow

**Use Case:** End-to-end audio for a promotional video

```python
def create_promo_video_audio():
    """Complete audio production for promotional video"""

    print("=== Promotional Video Audio Production ===\n")

    # 1. Background music
    print("Step 1: Generating background music...")
    bg_music = client.music_generation.compose(
        prompt="""Upbeat corporate background music with piano, light strings,
        and modern electronic elements. Professional and inspiring, perfect for
        company promotional video.""",
        music_length_ms=60000  # 1 minute
    )
    bg_path = save_audio(bg_music, "promo_background.mp3")

    # 2. Voiceover narration
    print("\nStep 2: Generating voiceover...")
    script = """
    At TechCorp, we're revolutionizing the way businesses connect with
    their customers. Our innovative solutions combine cutting-edge
    technology with human-centered design. Join thousands of companies
    who trust TechCorp to power their digital transformation.
    """

    voiceover = client.text_to_speech.convert(
        text=script,
        voice_id="pNInz6obpgDQGcFmaJgB",  # Adam - professional
        model_id="eleven_multilingual_v2"
    )
    vo_path = save_audio(voiceover, "promo_voiceover.mp3")

    # 3. Sound effects
    print("\nStep 3: Generating sound effects...")
    sfx_list = [
        ("logo_reveal", "company logo reveal sound, powerful whoosh and impact", 2.0),
        ("transition_1", "smooth transition whoosh, modern and clean", 1.0),
        ("transition_2", "elegant transition sound, professional", 1.0),
        ("end_sting", "ending logo sting, memorable and impactful", 3.0),
    ]

    for name, description, duration in sfx_list:
        audio = client.text_to_sound_effects.convert(
            text=description,
            duration_seconds=duration
        )
        save_audio(audio, f"promo_{name}.mp3")

    print("\n=== Production Complete ===")
    print("Audio files generated:")
    print(f"  - Background music: {bg_path}")
    print(f"  - Voiceover: {vo_path}")
    print(f"  - Sound effects: {len(sfx_list)} files")
    print("\nNext: Import into video editing software and mix")

# Run the workflow
create_promo_video_audio()
```

## Helper Functions Library

```python
# Useful helper functions for working with ElevenLabs

def list_available_voices(filter_by=None):
    """List all available voices with optional filtering"""
    voices = client.voices.get_all()

    for voice in voices.voices:
        if filter_by:
            labels_str = str(voice.labels).lower()
            if filter_by.lower() not in labels_str:
                continue

        print(f"\nName: {voice.name}")
        print(f"ID: {voice.voice_id}")
        print(f"Labels: {voice.labels}")
        print(f"Category: {voice.category}")

def test_voice_sample(text, voice_id):
    """Quick test of a voice with sample text"""
    audio = client.text_to_speech.convert(
        text=text,
        voice_id=voice_id,
        model_id="eleven_flash_v2_5"  # Fast for testing
    )
    save_audio(audio, f"voice_test_{voice_id}.mp3")

def estimate_character_count(text):
    """Estimate characters for quota planning"""
    return len(text)

def split_long_text(text, max_chars=5000):
    """Split long text at sentence boundaries"""
    sentences = text.replace('! ', '!|').replace('? ', '?|').replace('. ', '.|').split('|')
    chunks = []
    current_chunk = ""

    for sentence in sentences:
        if len(current_chunk) + len(sentence) < max_chars:
            current_chunk += sentence + " "
        else:
            if current_chunk:
                chunks.append(current_chunk.strip())
            current_chunk = sentence + " "

    if current_chunk:
        chunks.append(current_chunk.strip())

    return chunks
```

## Best Practices Summary

1. **Always handle errors** with try/except blocks
2. **Use appropriate models** for your use case
3. **Test voices** before generating long content
4. **Organize output files** with clear naming
5. **Monitor API usage** to avoid quota issues
6. **Save successful prompts** for reuse
7. **Process in batches** with retry logic
8. **Use helper functions** to avoid repetition

## Resources

- **API Documentation**: https://elevenlabs.io/docs
- **Python SDK**: https://github.com/elevenlabs/elevenlabs-python
- **Voice Library**: https://elevenlabs.io/voice-library
- **Usage Dashboard**: https://elevenlabs.io/app/usage
