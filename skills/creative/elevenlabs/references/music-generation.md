# ElevenLabs Music Generation Reference

## Overview

ElevenLabs Music Generation API creates original, royalty-free music from text descriptions. Generate vocal and instrumental tracks across various genres, styles, and moods for video backgrounds, game soundtracks, podcasts, and creative projects.

**Important:** Music generation is available only to **paid ElevenLabs subscribers**.

## Capabilities

### Text-to-Music
- Generate complete music tracks from text prompts
- Vocal and instrumental compositions
- Multiple genres and styles
- Customizable track duration
- Royalty-free output

### Composition Plans
- Structured music blueprints (JSON schemas)
- Define global and section-level styles
- Control music structure and progression
- Fine-grained composition control

### Customization
- Genre specification
- Mood and atmosphere control
- Instrument selection
- Track length configuration (milliseconds)
- Language selection (for vocals)

## API Usage

### Basic Setup

```python
import os
from elevenlabs.client import ElevenLabs
from pathlib import Path

# Initialize client
client = ElevenLabs(api_key=os.environ.get("ELEVENLABS_API_KEY"))

# Note: Requires paid subscription
```

### Simple Music Generation

```python
# Generate music from text prompt
try:
    audio = client.music_generation.compose(
        prompt="""Upbeat indie pop song with acoustic guitar, light drums,
        and cheerful melody. Modern and energetic feel.""",
        music_length_ms=30000  # 30 seconds
    )

    # Save to file
    output_path = Path("background_music.mp3")
    with output_path.open("wb") as f:
        for chunk in audio:
            f.write(chunk)

    print(f"Music saved: {output_path}")

except Exception as e:
    if "paid" in str(e).lower() or "subscription" in str(e).lower():
        print("Error: Music generation requires a paid subscription")
    else:
        raise
```

### Using Composition Plans

```python
# Create composition plan first (more control)
composition_plan = client.music_generation.composition_plan.create(
    prompt="""Progressive house electronic track with build-up, energetic drop,
    and calm outro section""",
    music_length_ms=60000  # 60 seconds
)

# Generate music from plan
audio = client.music_generation.compose(
    composition_plan=composition_plan
)

# Save music
output_path = Path("edm_track.mp3")
with output_path.open("wb") as f:
    for chunk in audio:
        f.write(chunk)
```

## Parameters

### compose() Method

**prompt** (string)
- Text description of desired music
- Mutually exclusive with `composition_plan`
- Include genre, instruments, mood, style
- Examples shown below in Prompt Engineering section

**music_length_ms** (integer)
- Track duration in milliseconds
- Typical range: 10,000ms (10s) to 120,000ms (120s)
- Longer tracks take more processing time
- Example: `30000` for 30-second track

**composition_plan** (object)
- Structured music blueprint
- Mutually exclusive with `prompt`
- Created via `composition_plan.create()`
- Provides fine-grained control over structure

## Prompt Engineering

### Effective Prompt Structure

**Formula:** `[Genre/Style] + [Instruments] + [Mood/Atmosphere] + [Context/Use Case] + [Additional Details]`

### Example Prompts

**Basic:**
```python
prompt = "Upbeat pop music"
```

**Better:**
```python
prompt = "Upbeat pop music with piano and drums, cheerful and energetic"
```

**Best:**
```python
prompt = """Upbeat indie pop song with acoustic guitar, light drums, and
cheerful piano melody. Modern and energetic feel, perfect for lifestyle
video background. Instrumental only."""
```

### Prompt Components

**1. Genre/Style**
- Be specific about musical genre
- Examples: "indie pop", "cinematic orchestral", "lo-fi hip hop", "progressive house"

**2. Instruments**
- List primary instruments
- Examples: "acoustic guitar, piano, soft drums"
- Can specify instrument roles: "lead guitar, bass, percussion"

**3. Mood/Atmosphere**
- Emotional quality
- Examples: "cheerful", "melancholic", "dramatic", "peaceful", "energetic"

**4. Structure (optional)**
- Mention sections if important
- Examples: "with build-up and drop", "intro, main theme, outro"

**5. Context/Use Case**
- Where music will be used
- Examples: "perfect for video background", "game menu music", "podcast intro"

**6. Vocal vs Instrumental**
- Specify explicitly
- "instrumental only" or "with vocals"
- For vocals: can specify language

### Genre-Specific Prompts

**Cinematic/Orchestral:**
```python
prompts = {
    "epic": """Epic cinematic orchestral music with dramatic strings, powerful
    brass, and heroic theme. Grand and inspiring, perfect for movie trailer.""",

    "emotional": """Emotional piano and string composition, melancholic and
    beautiful, touching atmosphere, perfect for dramatic scene.""",

    "adventure": """Adventure orchestral music with playful woodwinds, energetic
    strings, and uplifting brass. Exciting and whimsical."""
}
```

**Electronic:**
```python
prompts = {
    "edm": """Progressive house track with energetic synths, powerful bass drops,
    and building tension. Club atmosphere, dance music.""",

    "ambient": """Ambient electronic soundscape with ethereal pads, subtle textures,
    and peaceful atmosphere. Meditative and calming.""",

    "synthwave": """Retro synthwave music with nostalgic 80s synths, driving bass,
    and neon atmosphere. Cyberpunk aesthetic."""
}
```

**Acoustic/Organic:**
```python
prompts = {
    "folk": """Acoustic folk song with fingerstyle guitar, soft harmonica, and
    warm atmosphere. Rustic and intimate, campfire feel.""",

    "jazz": """Smooth jazz with piano, double bass, and gentle brushed drums.
    Sophisticated and relaxing, late night atmosphere.""",

    "classical": """Classical piano composition in romantic style, expressive and
    lyrical, flowing melody with rich harmonies."""
}
```

**Hip Hop/Beats:**
```python
prompts = {
    "lo-fi": """Chill lo-fi hip hop beats with jazz piano samples, vinyl crackle,
    and mellow drums. Relaxing study music atmosphere, instrumental.""",

    "trap": """Dark trap beat with heavy 808 bass, sharp hi-hats, and atmospheric
    pads. Moody and aggressive, urban sound.""",

    "boom-bap": """Classic boom-bap hip hop beat with dusty drums, jazz sample,
    and scratching. Old school 90s style, head-nodding groove."""
}
```

**Game Music:**
```python
prompts = {
    "menu": """Fantasy RPG game menu music with harp, soft strings, and magical
    atmosphere. Medieval and mysterious, looping background music.""",

    "battle": """Intense battle music with fast strings, powerful percussion, and
    dramatic brass. Action-packed and energetic.""",

    "exploration": """Calm exploration music with ambient pads, light piano, and
    nature sounds. Peaceful and atmospheric."""
}
```

**Video Backgrounds:**
```python
prompts = {
    "corporate": """Professional corporate background music with piano, soft strings,
    and gentle percussion. Positive and uplifting, business presentation.""",

    "vlog": """Upbeat lifestyle vlog music with ukulele, hand claps, and whistling.
    Happy and energetic, modern social media feel.""",

    "tutorial": """Calm tutorial background music with soft piano and light
    electronic elements. Focused and non-distracting."""
}
```

## Content Policy

### Restricted Content

**NOT Allowed:**
- Artist names (e.g., "in the style of The Beatles")
- Band names
- Song titles
- Album names
- Trademarked works
- Copyrighted material references

**Allowed:**
- Genre descriptions (e.g., "rock", "jazz", "classical")
- Style characteristics (e.g., "with guitar solos", "upbeat tempo")
- Mood and atmosphere (e.g., "energetic", "melancholic")
- Instrument specifications
- General musical terms

### Error Handling

```python
# Handle copyright-related errors
try:
    audio = client.music_generation.compose(
        prompt=prompt,
        music_length_ms=30000
    )
except Exception as e:
    error_msg = str(e).lower()

    if "bad_prompt" in error_msg:
        print("Prompt contains restricted content")
        print("Avoid artist names, bands, or copyrighted material")

        # API may provide suggested alternative
        if hasattr(e, 'suggestion'):
            print(f"Suggested alternative: {e.suggestion}")

    elif "paid" in error_msg or "subscription" in error_msg:
        print("Music generation requires paid subscription")

    else:
        print(f"Error: {e}")
```

## Duration Guidelines

### Short Clips (10-30 seconds)
- Background loops
- Transition music
- Short intros/outros
- Social media content

```python
short_music = [
    ("intro", "Upbeat podcast intro music, 10 seconds", 10000),
    ("transition", "Smooth transition music", 15000),
    ("logo", "Brand logo music sting", 5000)
]
```

### Medium Tracks (30-60 seconds)
- Video backgrounds
- Full intros/outros
- Commercial music
- Game loops

```python
medium_music = [
    ("video_bg", "Corporate presentation background", 45000),
    ("game_menu", "Game main menu music loop", 60000)
]
```

### Long Tracks (60-120 seconds)
- Complete songs
- Extended backgrounds
- Full compositions
- Podcast themes

```python
long_music = [
    ("full_song", "Complete indie pop song", 120000),
    ("ambient_track", "Long ambient soundscape", 180000)
]
```

## Use Cases

### Video Production

```python
video_music = {
    "opener": """Energetic opener music with electric guitar and drums,
    modern rock style, attention-grabbing intro.""",

    "background": """Subtle background music with soft piano and strings,
    calm and non-distracting, perfect for voice-over.""",

    "closer": """Uplifting outro music with full band, celebratory and
    positive, strong ending."""
}

for name, prompt in video_music.items():
    audio = client.music_generation.compose(
        prompt=prompt,
        music_length_ms=20000
    )
    save_audio(audio, f"video_{name}.mp3")
```

### Game Development

```python
game_tracks = [
    ("main_menu", """Fantasy RPG main menu music with orchestral strings,
    harp, and mystical atmosphere. Medieval and epic.""", 60000),

    ("combat", """Fast-paced battle music with intense drums, aggressive
    strings, and dramatic brass. Action-packed.""", 90000),

    ("peaceful", """Calm village music with flute, acoustic guitar, and
    gentle percussion. Peaceful and relaxing.""", 120000),

    ("victory", """Triumphant victory fanfare with brass, percussion, and
    uplifting melody. Celebratory.""", 15000)
]

for name, prompt, duration in game_tracks:
    audio = client.music_generation.compose(
        prompt=prompt,
        music_length_ms=duration
    )
    save_audio(audio, f"game_{name}.mp3")
```

### Podcast Production

```python
# Podcast theme and segments
podcast_music = {
    "theme": ("""Podcast theme music with modern electronic beats, synth
    melody, and professional sound. Catchy and memorable.""", 30000),

    "interlude": ("""Short interlude music with piano and soft pads,
    smooth transition between segments.""", 10000),

    "ad_break": ("""Upbeat ad break music with ukulele and hand claps,
    friendly and positive atmosphere.""", 8000),

    "outro": ("""Podcast outro music with fading electronic elements,
    calm ending theme.""", 15000)
}

for name, (prompt, duration) in podcast_music.items():
    audio = client.music_generation.compose(
        prompt=prompt,
        music_length_ms=duration
    )
    save_audio(audio, f"podcast_{name}.mp3")
```

### Content Creation

```python
# Social media and content backgrounds
content_music = [
    ("upbeat_vlog", """Happy vlog background music with ukulele, whistling,
    and hand claps. Fun and energetic, modern social media.""", 45000),

    ("tutorial_bg", """Calm tutorial background music with soft piano and
    light ambient pads. Non-distracting and focused.""", 120000),

    ("timelapse", """Inspiring timelapse music with building energy, piano
    and strings crescendo. Motivational.""", 60000)
]

for name, prompt, duration in content_music:
    audio = client.music_generation.compose(
        prompt=prompt,
        music_length_ms=duration
    )
    save_audio(audio, f"content_{name}.mp3")
```

## Composition Plans

### Understanding Composition Plans

Composition plans provide structured control over music generation:

**Global Styles:**
- Apply to entire composition
- Positive styles (what to include)
- Negative styles (what to avoid)

**Sections:**
- Individual segments of the track
- Section-specific styles
- Duration control per section
- Structural organization

### Creating Composition Plans

```python
# Generate composition plan from prompt
plan = client.music_generation.composition_plan.create(
    prompt="""Electronic dance music with intro build-up, energetic drop
    section, breakdown, and outro""",
    music_length_ms=90000
)

# Inspect the plan
print("Global Styles:", plan.global_styles)
print("Number of Sections:", len(plan.sections))
for i, section in enumerate(plan.sections):
    print(f"Section {i+1}: {section.duration_ms}ms")
```

### Using Composition Plans

```python
# Method 1: Generate plan first, then compose
plan = client.music_generation.composition_plan.create(
    prompt="Progressive house track with clear sections",
    music_length_ms=60000
)

audio = client.music_generation.compose(
    composition_plan=plan
)

# Method 2: Direct composition with prompt
audio = client.music_generation.compose(
    prompt="Progressive house track",
    music_length_ms=60000
)
# API generates plan automatically
```

## Best Practices

### 1. Clear Prompts
- Be specific about genre and style
- Mention key instruments
- Describe mood and atmosphere
- Specify vocal vs instrumental

### 2. Appropriate Duration
- Match duration to use case
- Shorter for loops and backgrounds
- Longer for complete compositions
- Consider processing time

### 3. Avoid Copyright
- Use generic genre descriptions
- Don't mention specific artists
- Focus on characteristics, not examples
- Use musical terminology

### 4. Iteration
- Generate multiple variations
- Try different prompt phrasings
- Adjust duration as needed
- Test in target application

### 5. Organization
- Name files descriptively
- Track successful prompts
- Organize by project/use case
- Document what worked

### 6. Rights & Usage
- Generated music is royalty-free
- Can be used commercially (with paid plan)
- Check terms of service for specifics
- Attribution not required but appreciated

## Advanced Techniques

### Generating Variations

```python
# Create variations of a theme
base_prompt = "Upbeat corporate background music with piano and strings"

variations = [
    f"{base_prompt}, energetic and fast-paced",
    f"{base_prompt}, calm and gentle",
    f"{base_prompt}, dramatic and powerful"
]

for i, prompt in enumerate(variations):
    audio = client.music_generation.compose(
        prompt=prompt,
        music_length_ms=30000
    )
    save_audio(audio, f"corporate_var{i+1}.mp3")
```

### Multi-Section Tracks

```python
# Use composition plan for structured multi-section music
complex_prompt = """
Cinematic trailer music with:
- Quiet mysterious intro with strings (10 seconds)
- Building tension with added percussion (15 seconds)
- Powerful climax with full orchestra (10 seconds)
- Dramatic ending (5 seconds)
"""

plan = client.music_generation.composition_plan.create(
    prompt=complex_prompt,
    music_length_ms=40000
)

audio = client.music_generation.compose(composition_plan=plan)
save_audio(audio, "trailer_music.mp3")
```

### Batch Generation

```python
def generate_music_library(tracks):
    """Generate multiple music tracks"""
    for name, prompt, duration in tracks:
        print(f"Generating: {name}...")

        try:
            audio = client.music_generation.compose(
                prompt=prompt,
                music_length_ms=duration
            )

            output_path = Path(f"music_{name}.mp3")
            with output_path.open("wb") as f:
                for chunk in audio:
                    f.write(chunk)

            print(f"  ✓ Saved: {output_path}")

        except Exception as e:
            print(f"  ✗ Error: {e}")

# Generate library
music_library = [
    ("upbeat", "Happy upbeat music", 30000),
    ("calm", "Calm peaceful music", 45000),
    ("dramatic", "Dramatic intense music", 60000)
]

generate_music_library(music_library)
```

## Troubleshooting

### Issue: "Paid subscription required"

**Solution:**
- Music generation is paid-tier only
- Upgrade account at https://elevenlabs.io/pricing
- Text-to-speech and sound effects available on free tier

### Issue: "Bad prompt" error

**Solution:**
- Remove artist names, band names, song titles
- Use generic genre and style descriptions
- Check API response for suggested alternative
- Rephrase using musical characteristics

### Issue: Generated music doesn't match description

**Solution:**
- Make prompt more specific
- Add more detail about instruments
- Specify mood and atmosphere clearly
- Try composition plan for more control

### Issue: Processing takes long time

**Solution:**
- Longer durations take more time
- Start with shorter tracks for testing
- Be patient (can take 30-60 seconds)
- Check API status if unusually slow

## API Limits

### Account Requirements
- **Free Tier:** Music generation not available
- **Paid Tiers:** Access to music generation API

### Rate Limits
- Subject to account tier limits
- Processing time increases with duration
- Monitor usage via dashboard

### Duration Limits
- Minimum: typically 5-10 seconds
- Maximum: varies by tier (usually up to 180 seconds)
- Longer tracks may require higher tiers

## Resources

- **Music Generation Docs**: https://elevenlabs.io/docs/cookbooks/music/quickstart
- **API Reference**: https://elevenlabs.io/docs/api-reference/music-generation
- **Pricing & Plans**: https://elevenlabs.io/pricing
- **Usage Dashboard**: https://elevenlabs.io/app/usage
- **Python SDK**: https://github.com/elevenlabs/elevenlabs-python
- **Community Showcase**: https://elevenlabs.io/showcase
